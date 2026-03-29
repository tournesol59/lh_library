#cp Frederic Peugny, LTS
import packcolsystwo as colsystwo

NSUB=1
NX=4
X=[-1.0, -0.66, -0.25, 0.25, 0.66, 1.0]
ALEFT=-1.0
ARIGHT=1.0
PARAM=[NX, ALEFT, ARIGHT, 'norm', NSUB]
colpts=colsystwo.ColsysBPoints(PARAM,X)

