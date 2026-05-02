import math as math
# generates the file code
numint=2;
intpoints=20;

with open('TheFileE.txt', 'w') as writer:
   #C1=t    C2=ut   C3=yt   C4=cd2   C5=cd1   C6=cd0
   for i in range(0,numint*intpoints+1):
      C1=i*0.03
      C2=0.
      C3=math.sin(2*C1)
      C4=1
      C5=0.5
      C6=16
# dont work     writer.write("{C1} \t {C2} \t {C3} \t {C4} \t {C5} \t {C6} \n ") 
      writer.write(f"{C1:02} \t {C2:02} \t {C3:02} \t {C4:02} \t {C5:02} \t {C6:02} \n ")
# writer close automatically inside the with clause

