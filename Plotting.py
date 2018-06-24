import matplotlib.pyplot as plt

class plotting:

    #def __init__(self,fp,fs,l_bl,t_bl,delta_est):


    @staticmethod
    def laminar(fp,fs,bl,delta_est,solutions,string,n):

        colors  = ["firebrick","sienna","black","lightseagreen","olivedrab","mediumvioletred"]

        fig = plt.figure(num=n, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
        if string == "vw":
            fig.suptitle("Blasius Laminar Flat plate Solution for different suction velocities (Re = %r)" %round(fs.Re,0), fontsize=18,fontweight="bold")
        elif string == "beta":
            fig.suptitle("Blasius Laminar Flat plate Solution for different pressure gradients/Falkner-Skan wedge angles (Re = %r)" %round(fs.Re,0), fontsize=18,fontweight="bold")
        
        # ## f-eta plot
        # ax = fig.add_subplot(231)
        # ax.plot(fp.eta,bl[str(0)].f ,'r-',linewidth=2,label='f')
        # ax.plot(fp.eta,bl[str(0)].fd ,'b-',linewidth=2,label='f´')
        # ax.plot(fp.eta,bl[str(0)].fdd ,'g-',linewidth=2,label='f´´' )
        # ax.set_title("Blasius Flat plate functions" + string + " = %r" %0)
        # ax.set_xlabel(r' $\eta$')
        # ax.set_xlim(0,4)
        # ax.set_ylim(0,5)
        # ax.legend()
        # ax.grid()

        #U-eta
        ax = fig.add_subplot(221)
        for idx,vw in enumerate(solutions):
            ax.plot(bl[str(vw)].fd,fp.eta,color = colors[idx],linestyle='-',linewidth=2,label=string + "=" +str(round(vw,2)))
        ax.set_title('Blasius Flat plate non dimensional Velocity U')
        ax.set_ylabel(r' $\eta$')
        ax.set_xlabel(r' $U/u_e$')
        ax.set_xlim(0,1)
        ax.set_ylim(0,6)
        ax.legend(loc = 'upper left')
        ax.grid()

        #V/(nu/)-eta
        ax = fig.add_subplot(222)
        for idx,vw in enumerate(solutions):
            ax.plot(fp.eta*bl[str(vw)].fd - bl[str(vw)].f,fp.eta,color = colors[idx],linestyle='-',linewidth=2,label=string + "=" +str(round(vw,2)))
        ax.set_title('Blasius Flat plate non dimensional Velocity V')
        ax.set_ylabel(r' $\eta$')
        ax.set_xlabel(r' $V/(\nu U/2x)^{0.5}$')
        if string == "vw":
            ax.set_xlim(-3,4)
        else:
            ax.set_xlim(0,2)
        ax.set_ylim(0,6)
        ax.legend(loc = 'upper left')
        ax.grid()

        ## boundary_layer thickness_plot
        ax = fig.add_subplot(223)

        for idx,vw in enumerate(solutions):    
            ax.plot(fp.x,bl[str(vw)].delta,color = colors[idx],linestyle='-',linewidth=2,label=r" $\delta$ " + string + " = "  +str(round(vw,2)))
        ax.set_title('Blasius Flat plate boundary layer')
        ax.set_xlabel(r' $x(m)$')
        ax.set_ylabel(r' $y(m)$')
        ax.set_xlim(0,1)
        ax.legend()
        ax.grid()

        # Cf Plot
        ax = fig.add_subplot(224)
        for idx,vw in enumerate(solutions): 
            ax.plot(fp.x, bl[str(vw)].cf,color=colors[idx],linestyle='-',linewidth=2,label= string + " = " +str(round(vw,2)) )
        ax.set_title('Blasius Flat plate Skin Friction')
        ax.set_xlabel(r' $x$')
        ax.set_ylabel(r' $C_f$')
        ax.set_xlim(0,1)
        ax.set_ylim(0,bl[str(0)].cf[0]/10)
        ax.legend(loc = 'upper right')
        ax.grid()

    @staticmethod
    def laminar_thermal(fp,fs,bl,delta_est,solutions,string,n):

        colors  = ["firebrick","sienna","black","lightseagreen","olivedrab","mediumvioletred"]

        fig = plt.figure(num=n, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')

        fig.suptitle("Blasius Laminar Flat plate Solution for different suction velocities (Re = %r)" %round(fs.Re,0), fontsize=18,fontweight="bold")
        
        # Theta-eta plot
        ax = fig.add_subplot(211)
        for idx,vw in enumerate(solutions):
            ax.plot(fp.eta,bl[str(vw)].theta,color = colors[idx],linestyle='-',linewidth=2,label= string + " = " +str(round(vw,2)))
        ax.set_title('Non dimensional Temperature for different Pr Numbers')
        ax.set_ylabel(r' $\theta$')
        ax.set_xlabel(r' $\eta$')
        ax.set_xlim(0,fp.eta[-1])
        ax.set_ylim(0,1)
        ax.legend(loc = 'upper right')
        ax.grid()

        # Nu-x Plot
        ax = fig.add_subplot(212)
        for idx,vw in enumerate(solutions): 
            ax.plot(fp.x, bl[str(vw)].Nu,color=colors[idx],linestyle='-',linewidth=2,label= string + " = " +str(round(vw,2)) )
        ax.set_title('Nusselt number Flat plate variation for different Pr numbers')
        ax.set_xlabel(r' $x$',fontsize=16)
        ax.set_ylabel(r' $Nu$',fontsize=16)
        ax.set_xlim(0,fp.x[-1])
        #ax.set_ylim(0,1)
        ax.legend(loc = 'upper left')
        ax.grid()

        # ax = fig.add_subplot(224)
        # ax.pcolor(fp.X,fp.Y, bl[str(0)].u)
        # ax.set_ylim(0,8*delta_est)
        # ax.plot(fp.x, bl[str(0)].delta,'b-',linewidth=2,label=r' $\delta$')
        # ax.plot(fp.x, bl[str(0)].delta_s,'r-',linewidth=2,label=r' $\delta^*$')
        # ax.plot(fp.x, bl[str(0)].delta_ss,'c-',linewidth=2,label=r' $\delta^{**}$')
        # ax.set_title("Blasius Flat plate boundary layer vw = %r" %0)
        # ax.set_xlabel(r' $x(m)$')
        # ax.set_ylabel(r' $y(m)$')
        # ax.set_xlim(0,1)
        # ax.legend()
        # ax.grid()

        fig.savefig('Images\\laminar_suction Re='+str(round(fs.Re,0))+'.png')

    @staticmethod
    def boundary_layer (fp,fs,bl,delta_est,ax,pcolor=False):

        ## Boundary Layer plot
        #fig3 = plt.figure(3)
        
        ax.plot(fp.x, bl.delta,'b-',linewidth=2,label=r' $\delta$')
        ax.plot(fp.x, bl.delta_s,'r-',linewidth=2,label=r' $\delta^*$')
        ax.plot(fp.x, bl.delta_ss,'c-',linewidth=2,label=r' $\delta^{**}$')
        #ax.title('Blasius Flat plate boundary layer')
        ax.xlabel(r' $x(m)$')
        ax.ylabel(r' $y(m)$')
        ax.xlim(0,1)
        ax.legend()
        ax.grid()



    @staticmethod
    def theta (fp,fs,bl):
        fig6 = plt.figure(6)
        plt.plot(fp.eta,bl.theta,linewidth=2,label="Pr = " + str(round(fs.Pr,3)))
        plt.title('Thermal Boundary Layer')
        plt.xlabel(r' $\eta$')
        plt.ylabel(r' $\theta$')
        plt.xlim(0,30)
        plt.ylim(0,1)
        plt.legend()
        plt.grid()

    @staticmethod 
    def drag_reynolds (Re,Cd_l,Cd_t_n,Cd_t_s,n):

        fig = plt.figure(num=n,figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
        plt.semilogx(Re, Cd_l,'k-',linewidth=2,label="Laminar")
        plt.semilogx(Re, Cd_t_n,'r-',linewidth=2,label="Turbulent")
        plt.semilogx(Re, Cd_t_s,'b-',linewidth=2,label="Transitional Model")
        plt.title('Drag Coefficient with Reynolds',fontsize=18)
        plt.xlabel(r' $Re$',fontsize=18)
        plt.ylabel(r' $C_d$',fontsize=18)
        plt.legend(prop={'size': 14})
        plt.grid()

        fig.savefig('Images\\Drag(Re).png')

    @staticmethod
    def xt_reynolds (Re,xt,n):

        fig = plt.figure(num=n,figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')
        plt.semilogx(Re, xt,'k-',linewidth=2)
        plt.title('Transition Point with Reynolds',fontsize=18)
        plt.xlabel(r' $Re$',fontsize=18)
        plt.ylabel(r' $x_t$',fontsize=18)
        plt.grid()

        fig.savefig('Images\\x_transition(Re).png')