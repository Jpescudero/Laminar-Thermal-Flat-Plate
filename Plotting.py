import matplotlib.pyplot as plt

class plotting:

    #def __init__(self,fp,fs,l_bl,t_bl,delta_est):


    @staticmethod
    def fblasius (fp,fs,bl):

        ## f-eta plot
        fig2 = plt.figure(2)
        plt.plot(fp.eta,bl.f ,'r-',linewidth=2,label='f')
        plt.plot(fp.eta,bl.fd ,'b-',linewidth=2,label='f´')
        plt.plot(fp.eta,bl.fdd ,'g-',linewidth=2,label='f´´' )
        plt.title('Blasius Flat plate functions')
        plt.xlabel(r' $\eta$')
        plt.xlim(0,4)
        plt.ylim(0,5)
        plt.legend()
        plt.grid()

    @staticmethod
    def boundary_layer (fp,fs,bl):
        ## Boundary Layer plot
        fig3 = plt.figure(3)
        #plt.pcolor(fp.X,fp.Y, bl.u)
        plt.plot(fp.x, bl.delta,'b-',linewidth=2,label=r' $\delta$')
        plt.plot(fp.x, bl.delta_s,'r-',linewidth=2,label=r' $\delta^*$')
        plt.plot(fp.x, bl.delta_ss,'c-',linewidth=2,label=r' $\delta^{**}$')
        plt.title('Blasius Flat plate boundary layer')
        plt.xlabel(r' $x(m)$')
        plt.ylabel(r' $y(m)$')
        plt.xlim(0,1)
        #plt.ylim(0,8*delta_est)
        #plt.colorbar()
        plt.legend()
        plt.grid()

    @staticmethod
    def cf (fp,fs,bl):
        ## Cf-x plot
        fig4 = plt.figure(4)
        plt.plot(fp.x, bl.cf,'k-',linewidth=2,label=r' $C_f$')
        plt.title('Blasius Flat plate Skin Friction')
        plt.xlabel(r' $x$')
        plt.ylabel(r' $C_f$')
        plt.xlim(0,1)
        #plt.ylim(0,1)
        plt.legend()
        plt.grid()

    @staticmethod
    def tw (fp,fs,bl):
        ## tw-x plot
        fig5 = plt.figure(5)
        plt.plot(fp.x, bl.tw,'k-',linewidth=2,label=r' $\tau_p$')
        plt.title('Blasius Flat plate wall shear stress at the wall')
        plt.xlabel(r' $x$')
        plt.ylabel(r' $\tau_p$')
        plt.xlim(0,1)
        plt.legend()
        plt.grid()

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
    def drag_reynolds (Re,Cd_l,Cd_t_n,Cd_t_s):
        fig7 = plt.figure(7)
        plt.semilogx(Re, Cd_l,'k-',linewidth=2,label="laminar")
        plt.semilogx(Re, Cd_t_n,'r-',linewidth=2,label="turbulent")
        plt.semilogx(Re, Cd_t_s,'b-',linewidth=2,label="Transitional Model")
        plt.title('Drag Coefficient with Reynolds')
        plt.xlabel(r' $Re$')
        plt.ylabel(r' $C_d$')
        plt.legend()
        plt.grid()

    @staticmethod
    def xt_reynolds (Re,xt):
        fig8 = plt.figure(8)
        plt.semilogx(Re, xt,'k-',linewidth=2)
        plt.title('Transition Point with Reynolds')
        plt.xlabel(r' $Re$')
        plt.ylabel(r' $x_t$')
        plt.grid()
