from PIL import Image
import numpy as np, sys
im=np.asarray(Image.open('sunstone_anagoria_3840.jpg').convert('L')).astype(float)
def unwrap(cx,cy,q,R0,R1,W):
    th=np.linspace(0,2*np.pi,W,endpoint=False); r=np.arange(R1,R0,-1.0)
    TT,RR=np.meshgrid(th,r)
    X=cx+RR*np.sin(TT); Y=cy-q*RR*np.cos(TT)
    xi=np.clip(X,0,im.shape[1]-2); yi=np.clip(Y,0,im.shape[0]-2)
    x0=xi.astype(int); y0=yi.astype(int); fx=xi-x0; fy=yi-y0
    return (im[y0,x0]*(1-fx)*(1-fy)+im[y0,x0+1]*fx*(1-fy)+im[y0+1,x0]*(1-fx)*fy+im[y0+1,x0+1]*fx*fy)
if __name__=='__main__':
    cx,cy,q=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])
    R0,R1,W=int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
    v=unwrap(cx,cy,q,R0,R1,W)
    out=Image.fromarray(np.clip(v,0,255).astype(np.uint8)); out.save(sys.argv[7])
    out.resize((2000,max(1,int(out.size[1]*2000/W)))).save(sys.argv[7].replace('.png','_small.png'))
