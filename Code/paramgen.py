import sys
import os

distance, beta2, N, mshapem, chirp0 = sys.argv[1::]
print distance, beta2, N, mshapem, chirp0
print str(sys.argv)
#name = "p_d"+str(distance)+"-b"+str(beta2)+"-N"+str(N)+"-sp"+str(mshapem)+"-ch"+str(chirp0)
param = open('paratemplate', 'r')
out =open("par-spec.py", 'w')
imp = ["import math\n","import numpy as np\n"]
out.writelines(imp)
text = ["distance="+str(distance)+"\n","beta2="+str(beta2)+"\n","N="+str(N)+"\n","mshape="+str(mshapem)+"\n","chirp0="+str(chirp0)+"\n"]
out.writelines(text)
out.write(param.read())
param.close
out.close
print distance, beta2, N, mshapem, chirp0
