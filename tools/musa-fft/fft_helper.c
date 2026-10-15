// wos-musa-fft — MUSA (Moore Threads) FFT helper daemon for Tetsuhiro-WOS.
//
// Why a helper process: the MUSA userspace driver (libmusart) reserves its
// device address window [4GB, ~873GB) inside the process at init time and
// aborts (musaErrorInvalidDevice, "invalid device ordinal") if any foreign
// VMA exists in that range. The Go runtime (linux/amd64) places its heap
// arena at 0xc000000000 (768GB) — inside that window — so NO Go process can
// initialize the device, while plain C processes can. This daemon runs in a
// plain C process (where the window is free) and serves FFTs to the Go side
// over a UNIX socket. See Tetsuhiro-WOS docs/notes-musa-rootcause.md.
//
// Protocol (little-endian, both directions; host = x86-64 on both sides):
//   request  header (24 bytes): u32 magic(0x46465557 "WUFF"), u32 op,
//                               u32 dir, u32 n, u64 bytes
//     op: 1 = 2-D FFT (n x n, row-major), 2 = 1-D FFT (length n)
//     dir: 0 = forward, 1 = inverse (unnormalized; the caller scales by 1/n
//          or 1/n^2, mirroring the pure-Go/WOS convention)
//     bytes = n*n*16 or n*16 (complex128 payload follows the header)
//   reply    header (16 bytes): u32 status (0=ok, else error), u32 rsvd,
//                               u64 bytes; then payload bytes when status==0.
//
// Build: see tools/musa-fft/build.sh (gcc -O2 -Wall, links libmufft/libmusart).
//
// NOTE (MUSA 3.1.0 pitfall): muFFT silently no-ops on host pointers — exec
// returns MUFFT_SUCCESS and leaves the buffers untouched. All transform
// buffers must live in device memory (musaMalloc + H2D/D2H), which is what
// both this daemon's serve path and --check do.
//
// Usage:
//   wos-musa-fft [--socket PATH] [--idle-exit SECONDS] [--verbose]
//   wos-musa-fft --check        # one-shot device self test, exit 0/1
//
// Clients are served concurrently: one detached thread per connection. Only
// the GPU section of a request (staging buffer, plan cache, device buffer,
// muFFT/MUSA calls) is serialized by g_gpu_mu, so a connection that stays
// open but idle never blocks other clients. The daemon exits after
// --idle-exit seconds with no active connection (0 = run until killed).
// The socket file is removed on exit.

#define _GNU_SOURCE
#include <mufft.h>
#include <musa_runtime.h>

#include <errno.h>
#include <poll.h>
#include <pthread.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>

#define MAGIC 0x46465557u /* "WUFF" */
#define OP_INFO 0u        /* info: reply payload = device name (no request payload) */
#define OP_FFT2D 1u
#define OP_FFT1D 2u
#define DIR_FWD 0u
#define DIR_INV 1u
#define MAX_BYTES_DEFAULT (16ull << 30) /* per-request payload cap (default): 16 GiB */

static const char *g_sock_path = "/tmp/wos-musa-fft.sock";
static int g_idle_exit = 600; /* seconds; 0 = never */
static int g_verbose = 0;
static volatile sig_atomic_t g_stop = 0;
static char g_dev_name[256] = "MUSA device";
static uint64_t g_max_bytes = MAX_BYTES_DEFAULT;

/* Concurrency: every client connection is served by its own detached thread.
   g_gpu_mu serializes only the GPU section of each request (staging buffer,
   plan cache, device buffer, muFFT/MUSA calls) — never the socket wait, so an
   idle client cannot starve the others. g_act_mu guards the live-client count
   and the last-activity timestamp used by the idle-exit check. */
static pthread_mutex_t g_gpu_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t g_act_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_active = 0;
static uint64_t g_last_ms = 0;

/* ------------------------------------------------------------------ */
/* muFFT plan cache                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    uint32_t op;
    uint32_t n;
    mufftHandle plan;
} plan_ent;

static plan_ent *g_plans = NULL;
static int g_nplans = 0, g_capplans = 0;

static mufftHandle get_plan(uint32_t op, uint32_t n) {
    for (int i = 0; i < g_nplans; i++) {
        if (g_plans[i].op == op && g_plans[i].n == n) return g_plans[i].plan;
    }
    mufftHandle h;
    mufftResult r;
    if (op == OP_FFT2D) {
        r = mufftPlan2d(&h, (int)n, (int)n, MUFFT_Z2Z);
    } else {
        r = mufftPlan1d(&h, (int)n, MUFFT_Z2Z, 1);
    }
    if (r != MUFFT_SUCCESS) {
        fprintf(stderr, "wos-musa-fft: mufft plan op=%u n=%u failed: %d\n", op, n, (int)r);
        return NULL;
    }
    if (g_nplans == g_capplans) {
        g_capplans = g_capplans ? g_capplans * 2 : 8;
        g_plans = realloc(g_plans, sizeof(plan_ent) * (size_t)g_capplans);
        if (!g_plans) {
            fprintf(stderr, "wos-musa-fft: plan cache OOM\n");
            exit(3);
        }
    }
    g_plans[g_nplans].op = op;
    g_plans[g_nplans].n = n;
    g_plans[g_nplans].plan = h;
    g_nplans++;
    return h;
}

/* Grow-only device buffer (requests are serialized; one buffer suffices). */
static void *g_dev = NULL;
static size_t g_dev_cap = 0;

/* Grow-only host staging buffer (see serve_one). */
static double *g_stage = NULL;
static size_t g_stage_cap = 0;

static void *dev_buf(size_t nb) {
    if (g_dev && nb <= g_dev_cap) return g_dev;
    if (g_dev) {
        musaFree(g_dev);
        g_dev = NULL;
        g_dev_cap = 0;
    }
    void *p = NULL;
    musaError_t e = musaMalloc(&p, nb);
    if (e != musaSuccess) {
        fprintf(stderr, "wos-musa-fft: musaMalloc(%zu) failed: %d\n", nb, (int)e);
        return NULL;
    }
    g_dev = p;
    g_dev_cap = nb;
    return g_dev;
}

/* ------------------------------------------------------------------ */
/* IO helpers                                                          */
/* ------------------------------------------------------------------ */

static ssize_t read_full(int fd, void *buf, size_t n) {
    char *p = buf;
    size_t got = 0;
    while (got < n) {
        ssize_t r = read(fd, p + got, n - got);
        if (r < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (r == 0) return (ssize_t)got; /* EOF */
        got += (size_t)r;
    }
    return (ssize_t)got;
}

static int write_full(int fd, const void *buf, size_t n) {
    const char *p = buf;
    size_t sent = 0;
    while (sent < n) {
        ssize_t r = write(fd, p + sent, n - sent);
        if (r < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        sent += (size_t)r;
    }
    return 0;
}

static uint64_t now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000u + (uint64_t)ts.tv_nsec / 1000000u;
}

/* ------------------------------------------------------------------ */
/* Request handling                                                    */
/* ------------------------------------------------------------------ */

static int serve_gpu(int fd, uint32_t op, uint32_t dir, uint32_t n, uint64_t bytes, uint64_t need);

static int serve_one(int fd) {
    uint8_t req[24];
    if (read_full(fd, req, sizeof(req)) != (ssize_t)sizeof(req)) return -1;

    uint32_t magic, op, dir, n;
    uint64_t bytes;
    memcpy(&magic, req + 0, 4);
    memcpy(&op, req + 4, 4);
    memcpy(&dir, req + 8, 4);
    memcpy(&n, req + 12, 4);
    memcpy(&bytes, req + 16, 8);

    uint8_t rep[16];
    memset(rep, 0, sizeof(rep));

    if (magic == MAGIC && op == OP_INFO && n == 0 && bytes == 0) {
        uint64_t nb = strlen(g_dev_name);
        memcpy(rep + 8, &nb, 8);
        if (write_full(fd, rep, sizeof(rep)) != 0 ||
            write_full(fd, g_dev_name, nb) != 0) {
            return -1;
        }
        return 0;
    }

    if (magic != MAGIC || (op != OP_FFT2D && op != OP_FFT1D) ||
        (dir != DIR_FWD && dir != DIR_INV) || n == 0 || bytes > g_max_bytes) {
        uint32_t st = 1;
        memcpy(rep, &st, 4);
        write_full(fd, rep, sizeof(rep));
        fprintf(stderr, "wos-musa-fft: bad request magic=%08x op=%u dir=%u n=%u bytes=%llu\n",
                magic, op, dir, n, (unsigned long long)bytes);
        return -1;
    }
    /* Expected payload size: 16 bytes per complex128 element. */
    uint64_t need = (op == OP_FFT2D) ? (uint64_t)n * (uint64_t)n * 16ull
                                     : (uint64_t)n * 16ull;
    if (bytes != need) {
        uint32_t st = 2;
        memcpy(rep, &st, 4);
        write_full(fd, rep, sizeof(rep));
        fprintf(stderr, "wos-musa-fft: size mismatch op=%u n=%u bytes=%llu need=%llu\n",
                op, n, (unsigned long long)bytes, (unsigned long long)need);
        return -1;
    }

    pthread_mutex_lock(&g_gpu_mu);
    int rc = serve_gpu(fd, op, dir, n, bytes, need);
    pthread_mutex_unlock(&g_gpu_mu);
    return rc;
}

/* GPU section of one request. The caller must hold g_gpu_mu. */
static int serve_gpu(int fd, uint32_t op, uint32_t dir, uint32_t n, uint64_t bytes, uint64_t need) {
    uint8_t rep[16];
    memset(rep, 0, sizeof(rep));
    (void)musaSetDevice(0);

    /* Grow-only staging buffer: reusing it avoids re-faulting up to
       hundreds of MB of fresh anonymous pages on every request. */
    if (bytes > g_stage_cap) {
        free(g_stage);
        g_stage = malloc(bytes ? bytes : 1);
        if (!g_stage) {
            g_stage_cap = 0;
            uint32_t st = 3;
            memcpy(rep, &st, 4);
            write_full(fd, rep, sizeof(rep));
            return -1;
        }
        g_stage_cap = bytes;
    }
    double *host = g_stage;
    if (read_full(fd, host, bytes) != (ssize_t)bytes) {
        return -1;
    }

    uint64_t t0 = now_ms();
    int ok = 1;
    mufftHandle plan = get_plan(op, n);
    void *dev = plan ? dev_buf(bytes) : NULL;
    if (plan && dev) {
        musaError_t e = musaMemcpy(dev, host, bytes, musaMemcpyHostToDevice);
        if (e != musaSuccess) {
            fprintf(stderr, "wos-musa-fft: H2D failed: %d\n", (int)e);
        } else {
            mufftResult r = mufftExecZ2Z(plan, (mufftDoubleComplex *)dev,
                                         (mufftDoubleComplex *)dev,
                                         dir == DIR_FWD ? MUFFT_FORWARD : MUFFT_INVERSE);
            if (r != MUFFT_SUCCESS) {
                fprintf(stderr, "wos-musa-fft: exec failed: %d\n", (int)r);
            } else {
                musaDeviceSynchronize();
                e = musaMemcpy(host, dev, bytes, musaMemcpyDeviceToHost);
                if (e != musaSuccess) {
                    fprintf(stderr, "wos-musa-fft: D2H failed: %d\n", (int)e);
                } else {
                    ok = 0;
                }
            }
        }
    }

    if (ok) {
        uint32_t st = 4;
        memcpy(rep, &st, 4);
        write_full(fd, rep, sizeof(rep));
        return -1; /* request failed; connection dropped, caller falls back */
    }

    uint64_t rbytes = bytes;
    memcpy(rep + 8, &rbytes, 8);
    if (write_full(fd, rep, sizeof(rep)) != 0 ||
        write_full(fd, host, bytes) != 0) {
        return -1;
    }
    if (g_verbose) {
        fprintf(stderr, "wos-musa-fft: op=%u dir=%u n=%u bytes=%llu %.1f ms\n",
                op, dir, n, (unsigned long long)bytes, (double)(now_ms() - t0));
    }
    (void)need;
    return 0;
}

static void handle_client(int cfd) {
    for (;;) {
        int rc = serve_one(cfd);
        if (rc != 0) break; /* EOF or error: drop client */
    }
    close(cfd);
}

static void touch_activity(void) {
    pthread_mutex_lock(&g_act_mu);
    g_last_ms = now_ms();
    pthread_mutex_unlock(&g_act_mu);
}

static void *client_thread(void *arg) {
    int cfd = (int)(intptr_t)arg;
    pthread_mutex_lock(&g_act_mu);
    g_active++;
    pthread_mutex_unlock(&g_act_mu);
    handle_client(cfd);
    pthread_mutex_lock(&g_act_mu);
    g_active--;
    g_last_ms = now_ms();
    pthread_mutex_unlock(&g_act_mu);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Self test                                                           */
/* ------------------------------------------------------------------ */

static int self_check(void) {
    int count = 0;
    musaError_t e = musaGetDeviceCount(&count);
    if (e != musaSuccess || count <= 0) {
        fprintf(stderr, "check: musaGetDeviceCount failed: err=%d count=%d\n", (int)e, count);
        return 1;
    }
    struct musaDeviceProp prop;
    e = musaGetDeviceProperties(&prop, 0);
    printf("device 0: %s (%zu MB)\n", e == musaSuccess ? prop.name : "?",
           e == musaSuccess ? prop.totalGlobalMem >> 20 : 0);

    const int n = 8;
    mufftHandle plan;
    if (mufftPlan1d(&plan, n, MUFFT_Z2Z, 1) != MUFFT_SUCCESS) {
        fprintf(stderr, "check: plan failed\n");
        return 1;
    }
    mufftDoubleComplex hbuf[8];
    for (int i = 0; i < n; i++) {
        hbuf[i].x = i;
        hbuf[i].y = 0;
    }
    /* Device memory is mandatory (host pointers are silently ignored by
       muFFT); this mirrors the serve path byte for byte. */
    mufftDoubleComplex *d = NULL;
    if (musaMalloc((void **)&d, sizeof(hbuf)) != musaSuccess) {
        fprintf(stderr, "check: musaMalloc failed\n");
        return 1;
    }
    if (musaMemcpy(d, hbuf, sizeof(hbuf), musaMemcpyHostToDevice) != musaSuccess) {
        fprintf(stderr, "check: H2D failed\n");
        return 1;
    }
    if (mufftExecZ2Z(plan, d, d, MUFFT_FORWARD) != MUFFT_SUCCESS) {
        fprintf(stderr, "check: forward failed\n");
        return 1;
    }
    musaDeviceSynchronize();
    if (musaMemcpy(hbuf, d, sizeof(hbuf), musaMemcpyDeviceToHost) != musaSuccess) {
        fprintf(stderr, "check: D2H failed\n");
        return 1;
    }
    /* DC term of x[i]=i must be sum(i) = n(n-1)/2 = 28. */
    double dc = hbuf[0].x;
    if (mufftExecZ2Z(plan, d, d, MUFFT_INVERSE) != MUFFT_SUCCESS) {
        fprintf(stderr, "check: inverse failed\n");
        return 1;
    }
    musaDeviceSynchronize();
    if (musaMemcpy(hbuf, d, sizeof(hbuf), musaMemcpyDeviceToHost) != musaSuccess) {
        fprintf(stderr, "check: D2H(2) failed\n");
        return 1;
    }
    musaFree(d);
    double worst = 0;
    for (int i = 0; i < n; i++) {
        double x = hbuf[i].x / n - i;
        double y = hbuf[i].y / n;
        double dd = x * x + y * y;
        if (dd > worst) worst = dd;
    }
    printf("forward DC term = %.6f (expect 28)\n", dc);
    printf("round-trip worst |err|^2 = %.3e\n", worst);
    int ok = (worst < 1e-24) && (dc > 27.999 && dc < 28.001);
    printf("check: %s\n", ok ? "OK" : "SUSPECT");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void on_sig(int sig) {
    (void)sig;
    g_stop = 1;
}

int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--socket") && i + 1 < argc) {
            g_sock_path = argv[++i];
        } else if (!strcmp(argv[i], "--idle-exit") && i + 1 < argc) {
            g_idle_exit = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--max-bytes") && i + 1 < argc) {
            g_max_bytes = strtoull(argv[++i], NULL, 0);
            if (g_max_bytes == 0) g_max_bytes = MAX_BYTES_DEFAULT;
        } else if (!strcmp(argv[i], "--verbose")) {
            g_verbose = 1;
        } else if (!strcmp(argv[i], "--check")) {
            return self_check();
        } else {
            fprintf(stderr,
                    "usage: %s [--socket PATH] [--idle-exit SECONDS] [--max-bytes N] [--verbose] [--check]\n",
                    argv[0]);
            return 2;
        }
    }

    {
        int count = 0;
        musaError_t e = musaGetDeviceCount(&count);
        if (e != musaSuccess || count <= 0) {
            fprintf(stderr, "wos-musa-fft: no MUSA device (err=%d count=%d)\n", (int)e, count);
            return 2;
        }
    }
    struct musaDeviceProp prop;
    const char *name = "MUSA device";
    if (musaGetDeviceProperties(&prop, 0) == musaSuccess) name = prop.name;
    snprintf(g_dev_name, sizeof(g_dev_name), "%s", name);

    signal(SIGTERM, on_sig);
    signal(SIGINT, on_sig);

    unlink(g_sock_path);
    int lfd = socket(AF_UNIX, SOCK_STREAM, 0);
    if (lfd < 0) {
        perror("socket");
        return 1;
    }
    struct sockaddr_un addr;
    memset(&addr, 0, sizeof(addr));
    addr.sun_family = AF_UNIX;
    if (strlen(g_sock_path) >= sizeof(addr.sun_path)) {
        fprintf(stderr, "wos-musa-fft: socket path too long\n");
        return 1;
    }
    strncpy(addr.sun_path, g_sock_path, sizeof(addr.sun_path) - 1);
    if (bind(lfd, (struct sockaddr *)&addr, sizeof(addr)) != 0) {
        perror("bind");
        return 1;
    }
    if (listen(lfd, 64) != 0) {
        perror("listen");
        return 1;
    }
    fprintf(stderr, "wos-musa-fft: ready (device=%s, socket=%s)\n", name, g_sock_path);

    touch_activity();
    while (!g_stop) {
        struct pollfd pfd = {lfd, POLLIN, 0};
        int pr = poll(&pfd, 1, 1000);
        if (pr < 0) {
            if (errno == EINTR) continue;
            break;
        }
        if (pr == 0) {
            int idle = 0;
            pthread_mutex_lock(&g_act_mu);
            idle = g_active == 0 && g_idle_exit > 0 &&
                   now_ms() - g_last_ms > (uint64_t)g_idle_exit * 1000ull;
            pthread_mutex_unlock(&g_act_mu);
            if (idle) {
                fprintf(stderr, "wos-musa-fft: idle exit\n");
                break;
            }
            continue;
        }
        int cfd = accept(lfd, NULL, NULL);
        if (cfd < 0) {
            if (errno == EINTR) continue;
            if (errno == EAGAIN || errno == EWOULDBLOCK) continue;
            perror("accept");
            break;
        }
        touch_activity();
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
        pthread_t tid;
        if (pthread_create(&tid, &attr, client_thread, (void *)(intptr_t)cfd) != 0) {
            perror("pthread_create");
            close(cfd);
        }
        pthread_attr_destroy(&attr);
    }

    close(lfd);
    unlink(g_sock_path);
    if (g_dev) musaFree(g_dev);
    return 0;
}
