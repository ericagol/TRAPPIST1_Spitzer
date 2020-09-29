

using PyPlot
using JLD2
using Statistics

@load "../../../data/T1_libration_laplace.jld2" eloft1 eloft2 eloft3
cp = ["C0","C1","C2","C3","C4","C5","C6","C7"]
pname=["b","c","d","e","f","g","h"]
nskip0=[1,1,8,8,8,8,8,8]
clf()
figure(figsize=(6.3,6))
i1=1; i20=[210,210,8233,8233,8353,8383,8083]
for j=1:7; 
i2 = i20[j]
nskip = nskip0[j]
ecos1 = eloft1[i1:nskip:i2,3,j].*cos.(eloft1[i1:nskip:i2,4,j])
esin1 = eloft1[i1:nskip:i2,3,j].*sin.(eloft1[i1:nskip:i2,4,j])
ecos2 = eloft2[i1:nskip:i2,3,j].*cos.(eloft2[i1:nskip:i2,4,j])
esin2 = eloft2[i1:nskip:i2,3,j].*sin.(eloft2[i1:nskip:i2,4,j])
ecos3 = eloft3[i1:nskip:i2,3,j].*cos.(eloft3[i1:nskip:i2,4,j])
esin3 = eloft3[i1:nskip:i2,3,j].*sin.(eloft3[i1:nskip:i2,4,j])
#plot(ecos1 ,esin1,".")
#plot(ecos2 ,esin2,".")
#plot(ecos3 ,esin3,".")
plot(ecos1 .- mean(ecos1),esin1 .- mean(esin1),color=cp[j],label=pname[j],linewidth=2)
#plot(ecos2 .- mean(ecos2),esin2 .- mean(esin2),".",color=cp[j],alpha=0.5)
#plot(ecos3 .- mean(ecos3),esin3 .- mean(esin3),".",color=cp[j],alpha=0.5)
xlabel(L"e \cos{\omega} - \langle e \cos{\omega}\rangle",fontsize=15)
ylabel(L"e \sin{\omega} - \langle e \sin{\omega}\rangle",fontsize=15)
end
legend(fontsize=10)
tight_layout()
axis([-0.01,0.01,-0.01,0.01])
savefig("../T1_evector_forced.pdf",bbox_inches="tight")
#read(stdin,Char)

clf()
figure(figsize=(6.3,6))
for j=7:-1:1; 
i2 = i20[j]
nskip = nskip0[j]
ecos1 = eloft1[i1:nskip:i2,3,j].*cos.(eloft1[i1:nskip:i2,4,j])
esin1 = eloft1[i1:nskip:i2,3,j].*sin.(eloft1[i1:nskip:i2,4,j])
ecos2 = eloft2[i1:nskip:i2,3,j].*cos.(eloft2[i1:nskip:i2,4,j])
esin2 = eloft2[i1:nskip:i2,3,j].*sin.(eloft2[i1:nskip:i2,4,j])
ecos3 = eloft3[i1:nskip:i2,3,j].*cos.(eloft3[i1:nskip:i2,4,j])
esin3 = eloft3[i1:nskip:i2,3,j].*sin.(eloft3[i1:nskip:i2,4,j])
plot(ecos1 ,esin1,color=cp[j],label=pname[j],linewidth=2)
plot(ecos2 ,esin2,color=cp[j],linewidth=2)
plot(ecos3 ,esin3,color=cp[j],linewidth=2)
xlabel(L"e \cos{\omega}",fontsize=15)
ylabel(L"e \sin{\omega}",fontsize=15)
end
legend(fontsize=12)
tight_layout()
axis([-0.02,0.01,-0.015,0.015])
savefig("../T1_evector.pdf",bbox_inches="tight")
