# Draws images/U6/pleiades_members_gaia.png from the Lodieu et al. 2019 Gaia DR2 Pleiades table (run in this folder; needs numpy and matplotlib).
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
rows=[l.rstrip('\n').split('\t') for l in open('pleiades_lodieu2019_gaiadr2.tsv') if l.strip() and not l.startswith('#')]
hdr=rows[0];D=[r for r in rows[1:] if len(r)==len(hdr) and r[1].strip().replace('.','').isdigit()]
ra=np.array([float(r[1]) for r in D]);de=np.array([float(r[2]) for r in D]);G=np.array([float(r[3]) if r[3].strip() else 21 for r in D]);dc=np.array([float(r[9]) for r in D]);plx=np.array([float(r[4]) for r in D])
ra0,de0=56.8711,24.1051  # Alcyone
def gnom(a,d):
  a=np.radians(a);d=np.radians(d);a0=np.radians(ra0);d0=np.radians(de0)
  c=np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0)
  x=np.cos(d)*np.sin(a-a0)/c; y=(np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(a-a0))/c
  return -np.degrees(x),np.degrees(y)   # east to the left, as on the sky
x,y=gnom(ra,de);mem=dc<=11.6
dist_pc=1000/np.median(plx[mem]);tid_deg=np.degrees(np.arctan(11.6/dist_pc))
named={'Alcyone':(56.8711,24.1051),'Atlas':(57.2906,24.0534),'Electra':(56.2189,24.1133),'Maia':(56.4567,24.3678),'Merope':(56.5816,23.9484),'Taygeta':(56.3021,24.4672),'Pleione':(57.2968,24.1367),'Celaeno':(56.2009,24.2894),'Sterope':(56.4768,24.5545)}
def size(g): return np.clip(1.5*(13-g),0.3,None)**1.55+0.25
BG='#0b1020'
fig=plt.figure(figsize=(16,9.4),dpi=125,facecolor=BG);ax=fig.add_axes([0.05,0.08,0.56,0.82],facecolor=BG)
nm=~mem&(np.abs(x)<7.5)&(np.abs(y)<7.5)
ax.scatter(x[nm],y[nm],s=size(G[nm])*0.8,c='#6f7d9a',alpha=0.55,lw=0,label='Candidates beyond the tidal radius (%d in view)'%nm.sum())
o=np.argsort(-G[mem])
ax.scatter(x[mem][o],y[mem][o],s=size(G[mem][o]),c='#f4f1e6',lw=0,alpha=0.95,label='Members within the tidal radius (%d)'%mem.sum())
ax.add_patch(Circle((0,0),tid_deg,fill=False,ec='#e0b25a',lw=1.2,ls='--'))
ax.text(tid_deg*0.71+0.15,-tid_deg*0.71-0.35,'tidal radius 11.6 pc  (%.1f°)'%tid_deg,color='#e0b25a',fontsize=11)
ax.add_patch(Rectangle((-1.1,-0.75),2.2,1.5,fill=False,ec='#8fb0e6',lw=1))
ax.set_xlim(-7.5,7.5);ax.set_ylim(-7.5,7.5);ax.set_aspect('equal')
for s in ax.spines.values(): s.set_color('#3a4560')
ax.tick_params(colors='#9aa6bf',labelsize=10);ax.set_xlabel('degrees east (left) / west (right) of Alcyone',color='#9aa6bf',fontsize=11);ax.set_ylabel('degrees north / south of Alcyone',color='#9aa6bf',fontsize=11)
ax.set_xticks(range(-6,7,2));ax.set_xticklabels([str(-v) if v else '0' for v in range(-6,7,2)])
ax.grid(color='#1c2540',lw=0.6)
ax.legend(loc='lower left',facecolor='#141b30',edgecolor='#3a4560',labelcolor='#d8deea',fontsize=10,markerscale=1.5)
fig.text(0.05,0.955,'The Pleiades (M45): the whole cluster',color='#f4f1e6',fontsize=19,weight='bold')
fig.text(0.05,0.925,'Gaia DR2 members, Lodieu et al. 2019 (A&A 628, A66). Dot size = brightness (Gaia G). Distance %.0f pc (%.0f light years). North up, east left.'%(dist_pc,dist_pc*3.2616),color='#9aa6bf',fontsize=11)
# inset: core
ia=fig.add_axes([0.645,0.30,0.335,0.52],facecolor=BG)
c=(np.abs(x)<1.1)&(np.abs(y)<0.75)
cm=c&mem;cn=c&~mem
ia.scatter(x[cn],y[cn],s=size(G[cn])*1.7,c='#6f7d9a',lw=0,alpha=0.6)
ia.scatter(x[cm],y[cm],s=size(G[cm])*1.7,c='#f4f1e6',lw=0)
for n,(a,d) in named.items():
  px,py=gnom(np.array([a]),np.array([d]));off={'Pleione':(0.0,0.07),'Atlas':(0.0,-0.09),'Alcyone':(0,-0.09),'Merope':(0,-0.09),'Electra':(0,-0.09),'Celaeno':(0.0,0.06),'Maia':(0,-0.09),'Taygeta':(0,0.06),'Sterope':(0,0.06)}[n]
  ia.text(px[0]+off[0],py[0]+off[1],n,color='#ffd98a',fontsize=10,ha='center',va='center',weight='bold')
mx,my=gnom(np.array([56.5816]),np.array([23.9484]));ia.scatter(mx,my,s=140,facecolors='none',edgecolors='#ffd98a',lw=1.2)
ia.set_xlim(1.1,-1.1) if False else ia.set_xlim(-1.1,1.1);ia.set_ylim(-0.75,0.75);ia.set_aspect('equal');ia.set_xticks([]);ia.set_yticks([])
for s in ia.spines.values(): s.set_color('#8fb0e6')
ia.set_title('The bright core, 2.2° × 1.5° (blue box at left): the Seven Sisters with Atlas and Pleione',color='#cfe0ff',fontsize=10.5)
fig.text(0.645,0.255,'Grey dots are Gaia kinematic candidates just outside the tidal radius. Four of the\nbrightest stars (Atlas, Electra, Maia, Taygeta) fall there only because Gaia DR2 measures\nthe parallax of very bright stars less precisely; they are long-established members.\nMerope has no Gaia DR2 parallax, so it is not in the list; its place is shown by a ring.\n\nData: Lodieu, N. et al. 2019, A&A 628, A66, VizieR J/A+A/628/A66; ESA Gaia DR2.',color='#9aa6bf',fontsize=10,va='top')
fig.savefig('pleiades_members_gaia.png',facecolor=BG);print('tid',tid_deg,'dist',dist_pc,'mem',mem.sum(),'nonmem in view',nm.sum())
