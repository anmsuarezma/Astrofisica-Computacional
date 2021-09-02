// Omelyan Position-Extended-Forest-Ruth-Like
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double a=0.;
const double G=1;
const double c=1;
const double M=1;
const double rs=2*G*M/c;
//const double N=4;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;

// ------------------------------------ Clase Geodesica ---------------------------------------------
class Geodesica{
 private:
  //    -inf < t < inf   ;   r < inf   ;    0 < theta < pi     ;     0 <= phi < 2*pi
  double t, r, theta, phi;
  double Vt, Vr, Vtheta, Vphi;
  double Ft, Fr, Ftheta, Fphi;
 public:
  void Inicie(double t0, double r0, double theta0, double phi0,
	      double Vt0, double Vr0, double Vtheta0, double Vphi0);
  void CalculeFuerza();
  void MuevaR(double dl, double constante);
  void MuevaV(double dl, double constante);
  //void AgregueFuerza(double Vt, double Vr, double Vtheta, double Vphi);
  double Get_t(){return t;};
  double Get_r(){return r;};
  double Get_theta(){return theta;};
  double Get_phi(){return phi;};
  double GetVt(){return Vt;};
  double GetVr(){return Vr;};
  double GetVtheta(){return Vtheta;};
  double GetVphi(){return Vphi;};
  //friend class Schw;
};

void Geodesica::Inicie(double t0, double r0, double theta0, double phi0,
		  double Vt0, double Vr0, double Vtheta0, double Vphi0){
  t=t0;    r=r0;    theta=theta0;    phi=phi0;
  Vt=Vt0;  Vr=Vr0;  Vtheta0=Vtheta;  Vphi=Vphi0;
  Ft=0; Fr=0; Ftheta=0; Fphi=0;
}

void Geodesica::CalculeFuerza(){
  double Sigma=r*r+a*a*cos(theta)*cos(theta);
  double Delta=r*r-rs*r+a*a;
  double A=(r*r+a*a)*Sigma+rs*a*a*r*sin(theta)*sin(theta);
  // Gammas para t
  double Gamma001=rs*(r*r+a*a)*(r*r-a*a*cos(theta)*cos(theta))/(2*Sigma*Sigma*Delta);
  double Gamma002=(-1)*rs*a*a*r*sin(theta)*cos(theta)/(Sigma*Sigma);
  double Gamma023=rs*a*a*a*r*sin(theta)*sin(theta)*sin(theta)*cos(theta)/(c*Sigma*Sigma);
  double Gamma013=rs*a*sin(theta)*sin(theta)*(a*a*cos(theta)*cos(theta)*(a*a-r*r)-r*r*(a*a+3*r*r))/(2*c*Sigma*Sigma*Delta);
  // Gammas para r
  double Gamma100=c*c*rs*Delta*(r*r-a*a*cos(theta)*cos(theta))/(2*Sigma*Sigma*Sigma);
  double Gamma103=(-1)*c*Delta*rs*a*sin(theta)*sin(theta)*(r*r-a*a*cos(theta)*cos(theta))/(2*Sigma*Sigma*Sigma);
  double Gamma111=(2*r*a*a*sin(theta)*sin(theta)-rs*(r*r-a*a*cos(theta)*cos(theta)))/(2*Sigma*Delta);
  double Gamma112=(-1)*a*a*sin(theta)*cos(theta)/Sigma;
  double Gamma122=(-1)*r*Delta/Sigma;
  double Gamma133=(Delta*sin(theta)*sin(theta)/(2*Sigma*Sigma*Sigma))*(-2*r*Sigma*Sigma+rs*a*a*sin(theta)*sin(theta)*(r*r-a*a*cos(theta)*cos(theta)));
  // Gammas para theta
  double Gamma200=(-1)*c*c*rs*a*a*r*sin(theta)*cos(theta)/(Sigma*Sigma*Sigma);
  double Gamma203=c*rs*a*r*(r*r+a*a)*sin(theta)*cos(theta)/(Sigma*Sigma*Sigma);
  double Gamma211=a*a*sin(theta)*cos(theta)/(Sigma*Delta);
  double Gamma212=r/Sigma;
  double Gamma222=(-1)*a*a*sin(theta)*cos(theta)/Sigma;
  double Gamma233=(-1)*(sin(theta)*cos(theta)/(Sigma*Sigma*Sigma))*(A*Sigma+(r*r+a*a)*rs*a*a*r*sin(theta)*sin(theta));
  // Gammas para phi
  double Gamma301=c*rs*a*(r*r-a*a*cos(theta)*cos(theta))/(2*Sigma*Sigma*Delta);
  double Gamma302=(-1)*c*rs*a*r*cos(theta)/(Sigma*Sigma*sin(theta));
  double Gamma323=cos(theta)*(Sigma*Sigma+rs*a*a*r*sin(theta)*sin(theta))/(sin(theta)*Sigma*Sigma);
  double Gamma313=(2*r*Sigma*Sigma+rs*(a*a*a*a*sin(theta)*sin(theta)*cos(theta)*cos(theta)-r*r*(Sigma+r*r+a*a)))/(2*Sigma*Sigma*Delta);

  //Schwarzschild
  /*
  double Gamma100=(c*c*rs*(r-rs))/(2*r*r*r), Gamma001=rs/(2*r*(r-rs)), Gamma111=-rs/(2*r*(r-rs));
  double Gamma212=1/r, Gamma313=1/r, Gamma122=-(r-rs);
  double Gamma323=cos(theta)/sin(theta), Gamma133=-(r-rs)*sin(theta)*sin(theta), Gamma233=-sin(theta)*cos(theta);
  */
  //Fuerza t
  Ft=(-1)*(2*Gamma001*Vt*Vr+2*Gamma002*Vt*Vtheta+2*Gamma023*Vtheta*Vphi+2*Gamma013*Vr*Vphi);
  //Fuerza r
  Fr=(-1)*(Gamma100*Vt*Vt+Gamma111*Vr*Vr+Gamma122*Vtheta*Vtheta+Gamma133*Vphi*Vphi+2*Gamma103*Vt*Vphi+2*Gamma112*Vr*Vtheta);
  //Fuerza theta
  Ftheta=(-1)*(2*Gamma212*Vr*Vtheta+Gamma233*Vphi*Vphi+Gamma200*Vt*Vt+2*Gamma203*Vt*Vphi+Gamma211*Vr*Vr+Gamma222*Vtheta*Vtheta);
  //Fuerza phi
  Fphi=(-1)*(2*Gamma313*Vr*Vphi+2*Gamma323*Vtheta*Vphi+2*Gamma301*Vt*Vr+2*Gamma302*Vt*Vtheta);
  
}

void Geodesica::MuevaR(double dl, double constante){
  t+=Vt*(constante*dl); r+=Vr*(constante*dl);
  theta+=Vtheta*(constante*dl); phi+=Vphi*(constante*dl);
}
void Geodesica::MuevaV(double dl, double constante){
  Vt+=Ft*(constante*dl); r+=Fr*(constante*dl);
  Vtheta+=Ftheta*(constante*dl); Vphi+=Fphi*(constante*dl);
}

// ------------------------------------ Kerr -------------------------------------
int main(){
  // El l seria el parametro invariante
  double l=10, dl=-1e-4, lmax=-10;
  double b=1; // Parametro de impacto b=L/E
  Geodesica Foton;
  // Coordenadas iniciales
  double t0=-1, r0=3., theta0=M_PI/2 , phi0= M_PI/3;
  // Parametros necesarios para terminos de la metrica
  double Sigma=r0*r0+a*a*cos(theta0)*cos(theta0);
  double Delta=r0*r0-rs*r0+a*a;
  // Terminos necesarios de la metrica en condiciones iniciales
  double grr=Sigma/Delta;
  double gphiphi=(r0*r0+a*a+rs*r0*a*a*sin(theta0)*sin(theta0)/Sigma)*sin(theta0)*sin(theta0);
  double gtphi=2*rs*r0*a*sin(theta0)*sin(theta0)/Sigma;
  double gtt=-(1-rs*r0/Sigma);
  // Velocidades iniciales
  double Vtheta0=0;
  double Vt0=0.;//(-gphiphi-b*gtphi)/(gphiphi*gtt-gtphi*gtphi);
  double Vphi0=0.;//(b*gtt+b*gtphi)/(gphiphi*gtt-gtphi*gtphi);
  double Vr0=sqrt((-gtt*Vt0*Vt0-gphiphi*Vphi0*Vphi0-2*gtphi*Vt0*Vphi0)/grr);
  // Coordenadas preuso-cartesianas
  double x, y, z;
  
  
  Foton.Inicie(t0, r0, theta0, phi0, Vt0, Vr0, Vtheta0, Vphi0);

  for(; Foton.Get_r()>1.5; l+=dl){
    
    // dt = t-t0 ; dr = r-r0; dphi = phi2-phi1
    //cout<<Foton.Get_t()<<" "<<Foton.Get_r()<<" "<<Foton.Get_theta()<<" "<<Foton.Get_phi()<<endl;
    //cout<<Foton.GetVt()<<" "<<Foton.GetVr()<<" "<<Foton.GetVtheta()<<" "<<Foton.GetVphi()<<endl;
    
    x=sqrt((Foton.Get_r())*(Foton.Get_r())-a*a)*sin(Foton.Get_theta())*cos(Foton.Get_phi()-phi0);
    y=sqrt((Foton.Get_r())*(Foton.Get_r())-a*a)*sin(Foton.Get_theta())*sin(Foton.Get_phi()-phi0);
    z=(Foton.Get_r())*cos(Foton.Get_phi()-phi0);
    
    cout<<x<<" "<<y<<" "<<z<<endl;
    
    Foton.MuevaR(dl,Zeta);
    Foton.CalculeFuerza(); Foton.MuevaV(dl,(1-2*Lambda)/2);
    Foton.MuevaR(dl,Xi);
    Foton.CalculeFuerza(); Foton.MuevaV(dl,Lambda);
    Foton.MuevaR(dl,1-2*(Xi+Zeta));
    Foton.CalculeFuerza(); Foton.MuevaV(dl,Lambda);
    Foton.MuevaR(dl,Xi);
    Foton.CalculeFuerza(); Foton.MuevaV(dl,(1-2*Lambda)/2);
    Foton.MuevaR(dl,Zeta);
  }

  return 0;
}





