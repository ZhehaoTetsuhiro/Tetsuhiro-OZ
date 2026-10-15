package main

import (
	"log"
	"os"
	"strings"
	"time"

	"tetsuhiro-oz/optics"
)

// The helpers below are the WOS start-up behaviour (formerly cmd/wos/main.go),
// now shared by `oz serve`: backend selection, scripted-element definition
// directories with live reload, and the -gen-go exporter.

// announceBackend selects and logs the compute backend. -gpu asks for the GPU
// and reports when this build or machine cannot provide it; otherwise only the
// status is printed (WOS_GPU=1 in the environment still enables the GPU, see
// the optics package).
func announceBackend(wantGPU bool) {
	if !wantGPU {
		log.Printf("计算后端: %s", optics.GPUBackendInfo())
		return
	}
	if optics.SetGPU(true) {
		log.Printf("GPU 加速已启用: %s", optics.GPUBackendInfo())
	} else {
		log.Printf("-gpu 已指定但不可用（未以 -tags cuda/-tags musa 构建，或机器上没有对应设备）；继续使用 CPU")
	}
}

// setElementDirs installs the default definition directories followed by the
// comma-separated extra directories given with -elements.
func setElementDirs(extra string) {
	dirs := optics.DefaultElementDirs()
	for _, d := range strings.Split(extra, ",") {
		if d = strings.TrimSpace(d); d != "" {
			dirs = append(dirs, d)
		}
	}
	optics.SetElementDirs(dirs)
}

// reloadElements rescans the definition directories, logs the result and
// returns the report.
func reloadElements() optics.DefinitionReport {
	rep := optics.ReloadElementDefinitions(nil)
	logElementReport(rep)
	return rep
}

// runGenGo prints the native Go element for one definition file.
func runGenGo(path, typeName string) error {
	def, err := optics.LoadElementDefinitionFile(path)
	if err != nil {
		return err
	}
	if typeName == "" {
		typeName = def.Name
	}
	src, err := optics.GenerateElementGo(def, typeName)
	if err != nil {
		return err
	}
	_, err = os.Stdout.WriteString(src)
	return err
}

// logElementReport prints one line per loaded definition plus every error and
// note, so both the startup log and a reload tell the same story.
func logElementReport(rep optics.DefinitionReport) {
	log.Printf("元件定义: 扫描 %d 个目录，已加载 %d 个脚本元件", len(rep.Dirs), len(rep.Loaded))
	for _, d := range rep.Loaded {
		log.Printf("  加载 %s（%s · %s · %d 个参数）← %s", d.Name, d.Label, d.Behavior, d.Params, d.Source)
	}
	for _, n := range rep.Skipped {
		log.Printf("  跳过 %s: %s", n.Name, n.Reason)
	}
	for _, e := range rep.Errors {
		log.Printf("  错误 %s: %s", e.Source, e.Message)
	}
}

// watchElements reloads the definition directories whenever their contents
// change. A plain fingerprint poll keeps the kernel dependency-free; it
// compares the whole definition set, so a deleted or overwritten definition is
// noticed as well as an edited or new one.
func watchElements() {
	dirs := optics.ElementDirs()
	last := optics.DefinitionDirsSignature(dirs)
	for {
		time.Sleep(2 * time.Second)
		sig := optics.DefinitionDirsSignature(dirs)
		if sig == last {
			continue
		}
		last = sig
		log.Printf("元件定义有变化，重新加载…")
		logElementReport(optics.ReloadElementDefinitions(nil))
	}
}
