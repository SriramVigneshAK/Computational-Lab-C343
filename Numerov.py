#Numerov Method to Calculate Eiegnstate of NH3 Inversion Process(Double Well Potential)
import numpy as np
import matplotlib.pyplot as plt

 
E = 0.00555005
m =  1

x = np.linspace(-60, 60, num=1000)
s = x[2]-x[1]

V=[]
for i in range(0,len(x)):
  a = (0.0105567-3.5627*10**(-5)*(x[i]**(2))+3.12338*(10**(-8))*(x[i]**(4)))
  V.append(a)

def wavefn(E):
  Gn = []
  for i in range(0,len(x)):
    a = m*(2*(V[i]) - 2*E)
    Gn.append(a)

  psi = []
  psi.append(0)
  psi.append(0.0001)

  for i in range(2,len(x)):
    b = (2*psi[i-1]-psi[i-2]+((5*Gn[i-1]*psi[i-1]*(s**2))/6)+((Gn[i-2]*psi[i-2]*(s**2))/12))/(1-(Gn[i]*(s**2))/12)
    psi.append(b)
  return psi

def NODE():
  node = []
  for i in range(0,len(psi)-1):
    if psi[i]*psi[i+1]<0:
      node.append(1)
    else:
      node.append(0)
  return node

def itera():
if sum(node)==0:
  if psi[1000]<0.000000001:
    print("im here")
    if psi[1000]>-0.00000001:
      print("The value of E ", E)
    else:
      print("im ttthere")
      E = E + 0.005
      psi = wavefn(E)
  else:
    print("im heeeeere")
    E = E - 0.005
    psi = wavefn(E)
else:
  E = E + 0.005
  print("hello")
  psi = wavefn(E)
  
psi = wavefn(E)
node = NODE()


print(psi)
print("The number of nodes are", sum(node))
print("The energy is", E)



plt.plot(x, V, label = "V(x)")
plt.plot(x, psi, label = "psi(x)")
plt.xlabel('x')
plt.ylabel('psi')
plt.title('NH3 Eigenstates')
 
# show a legend on the plot
plt.legend()
 
# function to show the plot
plt.show()






