package main

// legacyPage is the original single-page OZ front end (geometry + coating).
// It stays reachable at /oz/ until the GUI tab merge lands; after that it is
// removed and /oz/ redirects to the unified GUI.
const legacyPage = `<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="theme-color" content="#101418">
<meta name="color-scheme" content="dark">
<title>Tetsuhiro-OZ</title>
<style>
 *{box-sizing:border-box}
 html{-webkit-text-size-adjust:100%}
 body{margin:0;padding:16px 14px 40px;background:#101418;color:#dde3ea;font:16px/1.55 system-ui,-apple-system,"Segoe UI","Noto Sans SC","PingFang SC","Microsoft YaHei",sans-serif}
 .wrap{max-width:880px;margin:0 auto}
 h1{display:flex;align-items:baseline;flex-wrap:wrap;gap:.5rem;font-size:1.3rem;font-weight:600;margin:.2rem 0 .5rem}
 h1 .ver{color:#7c8a99;font-size:.8rem;font-weight:400}
 .lede{color:#aab6c4;font-size:.95rem;margin:.1rem 0 1rem}
 .status{color:#7c8a99;font-size:.85rem}
 .card{background:#141a20;border:1px solid #232d38;border-radius:10px;padding:14px;margin:0 0 16px}
 .card h2{font-size:1rem;font-weight:600;margin:0 0 .8rem;color:#e6edf5}
 .endpoint{display:block;margin-top:2px;color:#7c8a99;font-size:.78rem;font-weight:400}
 code{background:#1c232b;padding:.1rem .35rem;border-radius:4px;font-size:.85em}
 .grid{display:grid;grid-template-columns:1fr;gap:10px}
 @media (min-width:560px){.grid{grid-template-columns:repeat(auto-fit,minmax(150px,1fr))}}
 .field label{display:block;margin:0 0 4px;color:#9fb0c0;font-size:.85rem}
 input,textarea,button{font-family:inherit;font-size:16px}
 input,textarea{width:100%;background:#161c22;color:#dde3ea;border:1px solid #2a3540;border-radius:8px;padding:.6rem .7rem;min-height:44px}
 textarea{min-height:7.5rem;resize:vertical;line-height:1.45}
 input:focus,textarea:focus,button:focus-visible{outline:2px solid #2f6feb;outline-offset:1px}
 .actions{display:flex;flex-wrap:wrap;gap:10px;margin:.8rem 0 .2rem}
 button{background:#2f6feb;color:#fff;border:0;border-radius:8px;min-height:44px;padding:.65rem 1.05rem;cursor:pointer;flex:1 1 auto;touch-action:manipulation;-webkit-tap-highlight-color:transparent}
 button.secondary{background:#26313c;color:#dde3ea;border:1px solid #2a3540}
 button:hover{filter:brightness(1.08)}
 button:disabled{opacity:.45;cursor:not-allowed}
 pre{background:#161c22;border:1px solid #2a3540;border-radius:8px;padding:.7rem;margin:.5rem 0 0;overflow:auto;max-height:24rem;white-space:pre-wrap;word-break:break-word;font:13px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
 canvas{display:block;width:100%;height:clamp(150px,38vw,300px);background:#0c0f13;border:1px solid #2a3540;border-radius:8px;margin:.9rem 0 .2rem}
 details.out{margin:.7rem 0 0}
 details.out summary{cursor:pointer;color:#8fb7ff;font-size:.9rem;padding:.35rem 0}
 .err{color:#ff9b9b}
 footer{color:#5b6a78;font-size:.8rem;text-align:center;margin-top:6px}
 @media (max-width:480px){
   body{padding:12px 10px 32px}
   h1{font-size:1.2rem}
   .card{padding:12px;border-radius:9px}
   .actions button{flex:1 1 100%}
 }
</style></head><body>
<div class="wrap">
<h1>Tetsuhiro-OZ <span class="ver">v` + version + `</span></h1>
<p class="lede">几何光学 · 镜头自主设计 · 镀膜自主设计 —— 本地服务。<code>GET /api/version</code> 查看状态。<br><span class="status" id="apiStatus">正在检查服务…</span></p>

<section class="card">
 <h2>① 消色差双胶合自主设计 <span class="endpoint">POST /api/lens/design · POST /api/lens/eval</span></h2>
 <div class="grid">
  <div class="field"><label for="eflMm">有效焦距 EFL（mm）</label>
   <input type="number" id="eflMm" value="100" min="1" step="1" inputmode="decimal"></div>
  <div class="field"><label for="semiAper">半口径（mm）</label>
   <input type="number" id="semiAper" value="10" min="0.1" step="0.5" inputmode="decimal"></div>
  <div class="field"><label for="fieldDeg">视场角（°）</label>
   <input type="number" id="fieldDeg" value="0.5" min="0" step="0.1" inputmode="decimal"></div>
  <div class="field"><label for="crownMat">冠玻璃 crown</label>
   <input type="text" id="crownMat" value="N-BK7" spellcheck="false" autocapitalize="off"></div>
  <div class="field"><label for="flintMat">火石玻璃 flint</label>
   <input type="text" id="flintMat" value="F2" spellcheck="false" autocapitalize="off"></div>
 </div>
 <div class="actions">
  <button id="designBtn" onclick="designLens()">设计</button>
  <button id="evalBtn" class="secondary" onclick="evalLens()" disabled>评估当前设计</button>
 </div>
 <details class="out" id="designDetails"><summary>设计结果 JSON（POST /api/lens/design）</summary>
  <pre id="designOut"></pre></details>
 <details class="out" id="evalOutDetails"><summary>评估结果 JSON（POST /api/lens/eval）</summary>
  <pre id="evalOut"></pre></details>
</section>

<section class="card">
 <h2>② 膜系设计 <span class="endpoint">POST /api/coating/design · POST /api/coating/spectrum</span></h2>
 <textarea id="coatReq" spellcheck="false" autocapitalize="off">{"template":"broadband-ar","from_nm":420,"to_nm":680,"incident":"Air","substrate":"BK7"}</textarea>
 <div class="actions">
  <button onclick="designCoating()">设计</button>
  <button onclick="spectrum()">设计并画光谱</button>
 </div>
 <canvas id="plot" aria-label="反射率光谱图"></canvas>
 <details class="out" id="coatDetails"><summary id="coatSummary">镀膜 JSON 结果</summary>
  <pre id="coatOut"></pre></details>
</section>

<footer>Tetsuhiro-OZ · 单页前端 · 全部接口见上（JSON 折叠区可复制结果）</footer>
</div>

<script>
const $=function(id){return document.getElementById(id);};
const LAMBDA_FROM=400, LAMBDA_TO=750;
let lastPrescription=null, lastRows=null;

async function api(url,body){
  const r=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});
  const t=await r.text();
  let d; try{d=JSON.parse(t);}catch(e){d={raw:t};}
  if(!r.ok){throw new Error((d&&d.error)?d.error:('HTTP '+r.status));}
  return d;
}
function setOut(id,v){
  const el=$(id); el.classList.remove('err');
  el.textContent=(typeof v==='string')?v:JSON.stringify(v,null,2);
}
function setErr(id,msg){const el=$(id); el.classList.add('err'); el.textContent='错误：'+msg;}
function numVal(id){const v=parseFloat($(id).value); return isNaN(v)?0:v;}
function openDetails(id){const d=$(id); if(d){d.open=true;}}

function lensRequest(){
  return {efl_mm:numVal('eflMm'),
          crown:$('crownMat').value.trim()||'N-BK7',
          flint:$('flintMat').value.trim()||'F2',
          semi_aperture_mm:numVal('semiAper'),
          field_deg:numVal('fieldDeg')};
}
async function designLens(){
  const b=$('designBtn'); b.disabled=true;
  try{
    const d=await api('/api/lens/design',lensRequest());
    lastPrescription=d.prescription||null;
    $('evalBtn').disabled=!lastPrescription;
    setOut('designOut',d); openDetails('designDetails');
  }catch(e){setErr('designOut',e.message); openDetails('designDetails');}
  finally{b.disabled=false;}
}
async function evalLens(){
  if(!lastPrescription){setErr('evalOut','请先完成一次镜头设计'); openDetails('evalOutDetails'); return;}
  const b=$('evalBtn'); b.disabled=true;
  try{
    const d=await api('/api/lens/eval',lastPrescription);
    setOut('evalOut',d); openDetails('evalOutDetails');
  }catch(e){setErr('evalOut',e.message); openDetails('evalOutDetails');}
  finally{b.disabled=false;}
}
async function designCoating(){
  let req;
  try{req=JSON.parse($('coatReq').value);}
  catch(e){setErr('coatOut','膜系请求 JSON 解析失败：'+e.message); return null;}
  try{
    const d=await api('/api/coating/design',req);
    setOut('coatOut',d);
    $('coatSummary').textContent='镀膜 JSON（设计结果）';
    return d;
  }catch(e){setErr('coatOut',e.message); return null;}
}
async function spectrum(){
  const d=await designCoating();
  if(!d||!d.stack)return;
  try{
    const s=await api('/api/coating/spectrum',{stack:d.stack,from_nm:LAMBDA_FROM,to_nm:LAMBDA_TO,step_nm:2,angle_deg:0,pol:'u'});
    lastRows=s.rows||[];
    $('coatSummary').textContent='镀膜 JSON（设计 + 光谱 '+lastRows.length+' 点）';
    drawPlot();
  }catch(e){setErr('coatOut','光谱计算失败：'+e.message);}
}
async function checkVersion(){
  try{
    const r=await fetch('/api/version');
    const d=await r.json();
    $('apiStatus').textContent='服务在线 · OZ '+d.version+' · '+d.go;
  }catch(e){$('apiStatus').textContent='服务状态未知（离线？）';}
}

/* --- 光谱图：宽度跟随容器，按 devicePixelRatio 设置后备像素 --- */
function fitCanvas(){
  const cv=$('plot');
  const dpr=window.devicePixelRatio||1;
  const w=Math.max(1,cv.clientWidth), h=Math.max(1,cv.clientHeight);
  cv.width=Math.round(w*dpr); cv.height=Math.round(h*dpr);
  const c=cv.getContext('2d');
  c.setTransform(dpr,0,0,dpr,0,0);
  return {c:c,w:w,h:h};
}
function drawPlot(){
  const f=fitCanvas(), c=f.c, w=f.w, h=f.h;
  c.clearRect(0,0,w,h);
  const x0=Math.max(34,Math.round(w*0.085)), x1=w-10, y0=22, y1=h-22;
  c.font='11px system-ui,sans-serif';
  c.strokeStyle='#2a3540'; c.lineWidth=1;
  c.beginPath(); c.moveTo(x0,y0); c.lineTo(x0,y1); c.lineTo(x1,y1); c.stroke();
  c.textAlign='right'; c.fillStyle='#7c8a99';
  for(const v of [0,0.5,1]){
    const y=y1-v*(y1-y0);
    c.strokeStyle=(v===0)?'#2a3540':'#1c232b';
    c.beginPath(); c.moveTo(x0,y); c.lineTo(x1,y); c.stroke();
    c.fillText(Math.round(v*100)+'%',x0-5,y+4);
  }
  c.textAlign='center';
  for(let lam=LAMBDA_FROM;lam<=LAMBDA_TO;lam+=50){
    const x=x0+(lam-LAMBDA_FROM)/(LAMBDA_TO-LAMBDA_FROM)*(x1-x0);
    c.strokeStyle='#2a3540';
    c.beginPath(); c.moveTo(x,y1); c.lineTo(x,y1+4); c.stroke();
    c.fillStyle='#7c8a99'; c.fillText(String(lam),x,y1+15);
  }
  c.textAlign='left'; c.fillStyle='#7c8a99';
  c.fillText('R(λ) %',x0,y0-9);
  c.textAlign='right';
  c.fillText('λ (nm) · 400-750 nm · 0°',x1,y0-9);
  if(!lastRows||!lastRows.length){
    c.textAlign='center'; c.fillStyle='#5b6a78';
    c.fillText('点「设计并画光谱」绘制反射率曲线',(x0+x1)/2,(y0+y1)/2);
    return;
  }
  c.strokeStyle='#5ea1ff'; c.lineWidth=2; c.beginPath();
  let first=true;
  for(const p of lastRows){
    const x=x0+(p.lambda_nm-LAMBDA_FROM)/(LAMBDA_TO-LAMBDA_FROM)*(x1-x0);
    const y=y1-p.R*(y1-y0);
    if(first){c.moveTo(x,y); first=false;} else {c.lineTo(x,y);}
  }
  c.stroke();
}
function refit(){drawPlot();}
window.addEventListener('resize',refit);
window.addEventListener('orientationchange',refit);
if(window.ResizeObserver){new ResizeObserver(refit).observe($('plot'));}
drawPlot();
checkVersion();
</script>
</body></html>
`
