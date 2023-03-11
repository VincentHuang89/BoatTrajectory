
N=1.2

plt.figure(dpi=300,figsize=(16, 10), constrained_layout=True)
phy = radians(60)
gamma = radians(20)
for a in range(1,12):
    X=[]
    Y=[]
    Y1=[]
    f=[]
    f1=[]
    for k in np.arange(0.4,30,0.1):
        f.append(F(k, N))
        f1.append(F1(k, N))
        X.append(x(a, k, N))
        Y.append((y(a, k / sin(phy - gamma), N)))
        Y1.append(-(y(a, k / sin(phy + gamma), N)))

    Attenuation0 = cos(phy - gamma)**2
    Attenuation1 = cos(phy + gamma)**2
    # Attenuation0 = abs(sin(2*(phy - gamma)))
    # Attenuation1 = abs(sin(2*(phy + gamma)))
    plt.plot(X, (Y), color='r', alpha=Attenuation0, linewidth=3)
    plt.plot(X, (Y1), color='r', alpha=Attenuation1, linewidth=3)

print(Attenuation0,Attenuation1)
plt.title('Froude Number: '+str(N))
plt.text(0,17,r'$\varphi$={:.1f}$\degree$, $\alpha$={:.1f}$\degree$'.format(degrees(phy), degrees(gamma)),fontsize=25)
plt.xlim(-1,25)
plt.ylim(-20,20)
plt.xlabel('Time', fontsize=15)
plt.ylabel("Channel", fontsize=15)