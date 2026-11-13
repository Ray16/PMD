#!/usr/bin/env python3
"""Build the all-design image overview and an offline interactive 3D gallery."""
import argparse
import html
import json
import shutil
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True)
    root=ap.parse_args().project_root.resolve();folder=root/'RFD4-proteina-designs/analysis/active_sites'
    records=json.loads((folder/'active_site_data.json').read_text())
    fontroot=Path('/usr/share/fonts/truetype/dejavu')
    font=lambda size,bold=False:ImageFont.truetype(str(fontroot/('DejaVuSans-Bold.ttf' if bold else 'DejaVuSans.ttf')),size)
    for kind,title in [('active_site','RFD4 active sites: both candidate proton-acceptor routes'),
                       ('proton_routes','RFD4 close-ups: Asp route and ATP gamma-oxygen route')]:
        canvas=Image.new('RGB',(3600,2160),'white');draw=ImageDraw.Draw(canvas)
        draw.text((60,25),title,font=font(48,True),fill='#18334d')
        draw.text((60,95),'Purple dashed: MVAP OH → Asp O   |   Orange dashed: MVAP OH → closest nonbridging ATP gamma O',font=font(27),fill='#44546a')
        draw.text((60,137),'MVAP carbon: green; ATP carbon: gold; Asp: purple; Lys: blue; Mg: teal. Same ligand alignment across designs.',font=font(25),fill='#44546a')
        for i,r in enumerate(records):
            x=(i%4)*900;y=205+(i//4)*940
            draw.rounded_rectangle((x+15,y,x+885,y+915),radius=15,outline='#d9e1e8',width=2)
            m=r['metrics'];base=r['Asp'];draw.text((x+35,y+15),f"{r['group']} · Asp{base['resi']} / Lys{r['Lys']['resi']}",font=font(30,True),fill='#18334d')
            im=Image.open(folder/'images'/f"{r['group']}_{kind}.png").convert('RGB')
            im.thumbnail((850,695));canvas.paste(im,(x+25+(850-im.width)//2,y+60+(695-im.height)//2))
            draw.text((x+35,y+765),f"Asp O: {m['Asp_O_to_MVAP_O16_A']:.3f} Å",font=font(26,True),fill='#7b3fb2')
            draw.text((x+35,y+804),f"ATP {r['gamma_acceptor']['name']}: {m['gamma_O_to_MVAP_O16_A']:.3f} Å",font=font(26,True),fill='#d55e00')
            ids=[v.split('_')[-1] for v in r['candidate_ids'].split(';')]
            draw.text((x+35,y+853),'Candidate IDs: '+', '.join(ids),font=font(20),fill='#44546a')
        draw.text((60,2100),'Heavy-atom contacts show geometric access; they do not identify the proton acceptor or establish a lower reaction barrier.',font=font(26),fill='#44546a')
        canvas.save(folder/f'all8_{kind}_overview.png')
    assets=folder/'assets';assets.mkdir(exist_ok=True)
    src=root/'saturation_mutagenesis/results/structural_atlas/assets'
    for name in ['3Dmol-min.js','3Dmol-min.js.LICENSE.txt']:shutil.copy2(src/name,assets/name)
    cards=[]
    for r in records:
        m=r['metrics'];cards.append(f'''<article><h2>{r['group']} — Asp{r['Asp']['resi']}</h2>
        <a href="images/{r['group']}_active_site.png"><img src="images/{r['group']}_active_site.png" alt="{r['group']} active site"></a>
        <p><span class="asp">Asp oxygen: {m['Asp_O_to_MVAP_O16_A']:.3f} Å</span> · <span class="gamma">ATP {r['gamma_acceptor']['name']}: {m['gamma_O_to_MVAP_O16_A']:.3f} Å</span></p>
        <p class="ids">{html.escape(r['candidate_ids'].replace(';',', '))}</p>
        <a href="images/{r['group']}_proton_routes.png">Close-up PNG</a> · <a href="images/{r['group']}_active_site.pse">PyMOL session</a> · <a href="structures/{r['group']}.pdb">Aligned PDB</a></article>''')
    options=''.join(f'<option value="{r["group"]}">{r["group"]}: {r["candidate_ids"].replace("candidate_", "").replace(";", ", ")}</option>' for r in records)
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>RFD4 active sites and proton acceptors</title><style>
    body{font:16px system-ui,sans-serif;color:#18334d;max-width:1500px;margin:auto;padding:28px;background:#f7f9fb}h1{font-size:32px;margin-bottom:8px}h2{font-size:20px}
    .intro{max-width:1100px;line-height:1.6}.asp{color:#7b3fb2}.gamma{color:#d55e00}a{color:#176398}button,select{font:inherit;padding:9px;margin:6px;border:1px solid #bbc9d4;border-radius:6px;background:white}
    .viewerwrap{display:grid;grid-template-columns:2fr 1fr;gap:20px;margin:25px 0}#viewer{height:600px;position:relative;background:white;border-radius:12px;border:1px solid #d9e1e8}aside{padding:20px;background:white;border-radius:12px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;text-align:left;border-bottom:1px solid #e4e9ee}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}article{background:white;padding:18px;border-radius:12px;border:1px solid #d9e1e8}article img{width:100%;height:auto}.ids{font-size:13px;color:#667584;line-height:1.5}.note{padding:14px;border-left:4px solid #cf8b30;background:#fff7e9;line-height:1.6}@media(max-width:800px){.viewerwrap,.grid{grid-template-columns:1fr}#viewer{height:450px}}
    </style><h1>RFD4 active sites: both candidate proton-acceptor routes</h1>
    <p class="intro">Eight distinct structures are represented by 37 repeated files. Purple marks the contact from MVAP C3 hydroxyl oxygen to the nearest Asp carboxylate oxygen. Orange marks its contact to the nearest nonbridging ATP gamma-phosphate oxygen. All three ATP gamma-oxygen distances are reported below.</p>
    <p class="note"><b>These are heavy-atom distances, not mechanism probabilities.</b> No hydroxyl hydrogen, validated protonation/solvation model, or matched activation barriers are available for the designs. The nearby Asp and Lys atoms carry source coordinate-conditioning annotations. The full catalytic network needs separate validation.</p>
    <p><a href="../proton_acceptor_distances.png">Two-distance plot</a> · <a href="../proton_acceptor_distances.csv">Full measurements CSV</a> · <a href="all8_active_site_overview.png">All 8 sites overview</a> · <a href="all8_proton_routes_overview.png">All 8 close-ups</a></p>
    <select id="choice">__OPTIONS__</select><button id="reset">Reset view</button><button id="close">Reactive center</button><label><input type="checkbox" id="labels"> Atom labels</label>
    <section class="viewerwrap"><div id="viewer"></div><aside><h2 id="gtitle"></h2><div id="measurements"></div><p class="ids" id="ids"></p><p class="ids">Drag to rotate; scroll to zoom. MVAP carbon is green, ATP carbon gold, Asp purple, Lys blue, Mg teal. Red atoms are oxygen, blue atoms nitrogen.</p></aside></section>
    <h2>All eight active sites</h2><div class="grid">__CARDS__</div>
    <script src="assets/3Dmol-min.js"></script><script>const sites=__DATA__;
    const viewer=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'white'});
    let current=sites[0];const coord=a=>({x:a.xyz[0],y:a.xyz[1],z:a.xyz[2]});
    function atomSel(a){return {chain:a.chain,resi:a.resi,atom:a.name};}
    function show(){current=sites.find(s=>s.group===document.getElementById('choice').value)||sites[0];let s=current,m=s.metrics;
    viewer.clear();viewer.addModel(s.pdb,'pdb');viewer.setStyle({chain:'H'},{cartoon:{color:'#c5ced6',opacity:0.2}});
    for(let resi of s.display_residues)viewer.setStyle({chain:'H',resi:resi},{stick:{radius:0.13,colorscheme:'Jmol'}});
    viewer.setStyle({chain:'L'},{stick:{radius:0.18,colorscheme:{prop:'elem',map:{C:'#269d62',O:'red',P:'#deb03d',N:'blue'}}}});
    viewer.setStyle({chain:'T'},{stick:{radius:0.17,colorscheme:{prop:'elem',map:{C:'#dbad39',O:'red',P:'#dbad39',N:'blue'}}}});
    viewer.setStyle({chain:'H',resi:s.Asp.resi},{stick:{radius:0.17,colorscheme:{prop:'elem',map:{C:'#7b3fb2',O:'red',N:'blue'}}}});
    viewer.setStyle({chain:'H',resi:s.Lys.resi},{stick:{radius:0.16,colorscheme:{prop:'elem',map:{C:'#4078ad',O:'red',N:'blue'}}}});
    viewer.setStyle({elem:'Mg'},{sphere:{radius:0.55,color:'#27aaa5'}});
    viewer.addLine({start:coord(s.hydroxyl),end:coord(s.Asp),color:'#7b3fb2',dashed:true,linewidth:3});
    viewer.addLine({start:coord(s.hydroxyl),end:coord(s.gamma_acceptor),color:'#d55e00',dashed:true,linewidth:3});
    if(document.getElementById('labels').checked){for(let [a,text] of [[s.hydroxyl,'MVAP O16'],[s.Asp,'Asp'+s.Asp.resi+' '+s.Asp.name],[s.gamma_acceptor,'ATP '+s.gamma_acceptor.name],[s.gamma_P,'ATP Pgamma']])viewer.addLabel(text,{position:coord(a),fontSize:13,fontColor:'#18334d',backgroundColor:'white',backgroundOpacity:0.8,borderThickness:0});}
    viewer.zoomTo({hetflag:true});viewer.render();document.getElementById('gtitle').textContent=s.group+' — candidate acceptor contacts';
    let rows=[['Asp '+s.Asp.resi+' '+s.Asp.name,m.Asp_O_to_MVAP_O16_A],['ATP IAT:O13',m.IAT_O13_to_MVAP_O16_A],['ATP IAT:O16',m.IAT_O16_to_MVAP_O16_A],['ATP IAT:O17',m.IAT_O17_to_MVAP_O16_A],['Glu oxygen (nearest)',m.Glu_O_to_MVAP_O16_A],['His nitrogen (nearest)',m.His_N_to_MVAP_O16_A],['Attack O16–Pgamma',m.attack_O_to_Pgamma_A]];
    document.getElementById('measurements').innerHTML='<table><tr><th>From MVAP O16 to</th><th>Distance</th></tr>'+rows.map(r=>'<tr><td>'+r[0]+'</td><td>'+r[1].toFixed(3)+' Å</td></tr>').join('')+'</table>';
    document.getElementById('ids').textContent=s.candidate_ids.replaceAll(';',', ');}
    document.getElementById('choice').onchange=show;document.getElementById('labels').onchange=show;document.getElementById('reset').onclick=show;
    document.getElementById('close').onclick=()=>{viewer.zoomTo({chain:'L',atom:'O16'});viewer.zoom(0.65);viewer.render();};show();
    </script></html>'''
    document=document.replace('__OPTIONS__',options).replace('__CARDS__',''.join(cards)).replace('__DATA__',json.dumps(records).replace('</','<\\/'))
    (folder/'index.html').write_text(document)
    print('Built both overview PNGs and the offline all-design 3D gallery.')


if __name__=='__main__':main()
