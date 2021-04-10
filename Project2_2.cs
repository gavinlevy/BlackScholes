using System;
class Project2{
    static void Main(){
        option();
        
    }
    
    static void option(){

        BSMOption p = new BSMOption(100, 1, 0.07,0.02,0.2);
        p.print();
        double s0=88;
        double t0=0;
        double Vtarget=18;
        double Vtol=0.01;
        p.Valuation(s0,t0);
        p.impliedVolatility(s0,t0,Vtarget,Vtol);
    }
}

public class BSMOption{
    
    private double SK=0;              
    private double T=0;             
    private double r=0;            
    private double q=0;          
    private double sigma=0;
    
    
    public BSMOption(double strike, double exTime, double rfRate,  double yield, double volatility) {
        SK = strike;
        T = exTime;
        r = rfRate;
        q = yield;
        sigma = volatility;
    }
    
    public double getK() { return SK; }
    public double getT() { return T; }
    public double getR() { return r; }
    public double getQ() { return q; }
    public double getSigma() { return sigma; }

    
    public void print(){
        Console.WriteLine($"strike price = {SK}");
        Console.WriteLine($"expiration time = {T}");
        Console.WriteLine($"risk-free rate = {r}");
        Console.WriteLine($"continuous dividend yield = {q}");
        Console.WriteLine($"volatility = {sigma}");

    }
    
    public void Valuation(double s0, double t0){
        double t = T-t0;
        double d1 = Math.Log10(s0/SK)+t*(r-q+(sigma*sigma/2));
        double d2 = d1-sigma*Math.Sqrt(t);
        double nd1 = RandLib.NormalCDF(d1);
        double nd2 = RandLib.NormalCDF(d2);
        double callFV=s0*Math.Exp(-q*t)*nd1-SK*Math.Exp(-r*t)*nd2;
        double putFV=SK*Math.Exp(-r*t)*RandLib.NormalCDF(-d2)-s0*Math.Exp(-q*t)*RandLib.NormalCDF(-d1);
        
        //delta
        double callDelta = Math.Exp(-q*t)*nd1;
        double putDelta = Math.Exp(-q*t)*(nd1-1);
        
        //gamma
        double gamma = (Math.Exp(-q*t)/(s0*sigma*Math.Sqrt(t)))*(1/Math.Sqrt(2*3.14))*Math.Exp(-d1*d1/2);
        
        //vega
        double vega = 0.01*s0*Math.Exp(-q*t)*Math.Sqrt(t)*(1/Math.Sqrt(2*Math.PI))*Math.Exp(-d1*d1/2);
        
        //Rho
        double callRho = 0.01*SK*t*Math.Exp(-r*t)*nd2;
        double putRho = -0.01*SK*t*Math.Exp(-r*t)*RandLib.NormalCDF(-d2);
        
        //theta
        double callTheta = 1/T*(-(s0*sigma*Math.Exp(-q*t)*Math.Exp(-d1*d1/2)/(2*Math.Sqrt(t)*Math.Sqrt(2*3.14)))-r*SK*Math.Exp(-r*t)*nd2+q*s0*Math.Exp(-q*t)*nd1);
        double putTheta = 1/T*(-(s0*sigma*Math.Exp(-q*t)*Math.Exp(-d1*d1/2)/(2*Math.Sqrt(t)*Math.Sqrt(2*3.14)))+r*SK*Math.Exp(-r*t)*RandLib.NormalCDF(-d2)-q*s0*Math.Exp(-q*t)*RandLib.NormalCDF(-d1));
        
        Console.WriteLine("call fair value ="+callFV);
        Console.WriteLine("put fair value ="+putFV);
        Console.WriteLine("call Delta ="+callDelta);
        Console.WriteLine("put Delta ="+putDelta);
        Console.WriteLine("gamma ="+gamma);
        Console.WriteLine("vega ="+vega);
        Console.WriteLine("call Rho ="+callRho);
        Console.WriteLine("put Rho ="+putRho);
        Console.WriteLine("call Theta ="+callTheta);
        Console.WriteLine("put Theta ="+putTheta);
        
        
        double pt=(callFV-putFV)-(s0*Math.Exp(-q*t)-SK*Math.Exp(-r*t));
        Console.WriteLine("put-call parity: "+pt);
   
    }
    public void impliedVolatility(double s0, double t0, double Vtarget, double Vtol){
        double t = T-t0;
        
        double[] x=new double[100];
        double[] si=new double[300];          

        si[0]=0.001;
        
        for(int i=0;i<100;i++){
            double d1 = Math.Log10(s0/SK)+t*(r-q+(si[i]*si[i]/2));
            double d2 = d1-si[i]*Math.Sqrt(t);
            double nd1 = RandLib.NormalCDF(d1);
            double nd2 = RandLib.NormalCDF(d2);
            double callFV = s0*Math.Exp(-q*t)*nd1-SK*Math.Exp(-r*t)*nd2;
            double putFV = SK*Math.Exp(-r*t)*RandLib.NormalCDF(-d2)-s0*Math.Exp(-q*t)*RandLib.NormalCDF(-d1);
            double vega = 0.01*s0*Math.Exp(-q*t)*Math.Sqrt(t)*(1/Math.Sqrt(2*Math.PI))*Math.Exp(-d1*d1/2);
            double gamma = (Math.Exp(-q*t)/(s0*sigma*Math.Sqrt(t)))*(1/Math.Sqrt(2*3.14))*Math.Exp(-d1*d1/2);

            x[i]=callFV;
            si[i+1]=si[i]+(-callFV+Vtarget)*gamma;

            
        }
        int count=0;
        for(int i=1;i<100;i++){
            count=count+1;
            if(si[i]==si[i+1]){
                break;
            }
            
        }
        
        Console.WriteLine("implied volatility="+si[count]*100+"%");
        
    }

}

sealed class RandLib
{
    private Random r;
    
    public static double NormalPDF(double x){
        return Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);
    }
    
    public static double NormalCDF(double x){
        if (x == 0.0) return 0.5;
        if (x < 0.0) return 1.0 - NormalCDF(-x);
        double[] q = { 0.6716419, 0.219386786, -0.53563782,
        1.6829937, -1.82639342, 1.7255361 };
        double t = 1.0 / (1.0 + q[0] * x);
        double sum = t * (q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * (q[5])))));
        double Phi = 1.0 - NormalPDF(x) * sum;
        return Phi;
    }

}
