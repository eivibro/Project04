from pylab import *
A = loadtxt('Program/build-Project4-Desktop-Debug/output.txt')
print A
energy = A[:,1]
temperature = A[:,0]
plot(temperature, energy)

J = 1
k = 1.381e-23
T = linspace(0.1, 20, 100)
E = -8*J*sinh(8/T)/(cosh(8/T)+3)
# Cv = 64/T*(1+3*cosh(8/T))/(cosh(8/T)+3)**2
# plot(T,Cv)
# ylabel('$C_V/k_B$', fontsize = 16)
# show()

plot(T,E,'r')
# ylabel('Energy')
show()
