
#ifndef PHYSICALMODEL_H_INCLUDED
#define PHYSICALMODEL_H_INCLUDED
#include <iostream>
#include <vector>
#include <numeric>
#include <random>

std::random_device rd;
std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
std::uniform_real_distribution<float> uni(-0.2,0.2); // guaranteed unbiased

std::vector<double> gl_RC;
std::vector<double> gl_spatial;
std::vector<double> gl_M1;
std::vector<double> gl_M2;
std::vector<double> gl_M3;

std::vector<double> gl_mem1;
std::vector<double> gl_mem2;
std::vector<double> gl_mem3;


////////////////////////////////////
// Mebrane class
class Membrane {
    float Lx,Ly,T,sigma0,sigma1_m,E,H,volDens;
    double* Am; //double* Bmf;
    
public:
    int Modes;
    float gamma;
    std::vector<double> spatial, M1,M2,M3,RC,mem1,mem2,mem3, MM1, MM2,MM3;
    Membrane(float _Lx,float _Ly,float _T,float _sigma0,float _sigma1_m,float _E,float _H,float _volDens,int _Modes,float _gamma) :
    Lx(_Lx),Ly(_Ly),T(_T),sigma0(_sigma0), sigma1_m(_sigma1_m),E(_E),H(_H),volDens(_volDens),Modes(_Modes),gamma(_gamma)
    {
        
        spatial.resize(Modes*Modes); gl_spatial.resize(Modes*Modes);
        M1.resize(Modes*Modes); M2.resize(Modes*Modes); M3.resize(Modes*Modes);
        gl_M1.resize(Modes*Modes); gl_M2.resize(Modes*Modes); gl_M3.resize(Modes*Modes);
        mem1.resize(Modes*Modes); mem2.resize(Modes*Modes); mem3.resize(Modes*Modes);
        gl_mem1.resize(Modes*Modes); gl_mem2.resize(Modes*Modes); gl_mem3.resize(Modes*Modes);
        MM1.resize(Modes*Modes); MM2.resize(Modes*Modes); MM3.resize(Modes*Modes);
        RC.resize(200); gl_RC.resize(200);
        
        Am  = (double*)calloc(Modes*Modes,sizeof(double));

        for (int i=0; i < 200;i++)
            RC[i] = 5 * (1-cos((float_Pi/200) * ( double(i) / 44100.0) )); // Raised cosine excitation signal
        

        
        float x_0 = 0.2; float y_0 = 0.7, nu = 0.3, D = (H*H*H*E) / ((12*(1-nu*nu))) , mu = (H * volDens);
        float k = 1/44100;
        float* m;
        float* n;
        m = (float*)calloc(Modes*Modes,sizeof(float)); n = (float*)calloc(Modes*Modes,sizeof(float));
        
        for (int i=0; i < Modes; i++){
            for (int j=0; j < Modes; j++){
                m[i*Modes+j] = i+1;
                n[i*Modes+j] = j+1;}}
        
        for (int j = 0; j < Modes*Modes; j++){
            mem1[j] = 0; mem2[j] = 0; mem3[j] = 0;
            Am[j]  = pow( (float_Pi * m[j]) / Lx ,2) + pow((float_Pi*n[j]) / Ly ,2) ;
            spatial[j] = sin(float_Pi*double(m[j])*(x_0/Lx)) * sin(float_Pi*double(n[j])*(y_0/Ly));
            // Current displacement
            M1[j] = (2.0 - double(5.1419e-10) * double( ( T * Am[j] + D * Am[j] * Am[j] ) / mu)) / (1 + double(k*(sigma0 + sigma1_m * Am[j])));
            // Previous displacement
            M2[j] = double((-1+ k*(sigma0 + sigma1_m * Am[j] ) ) / (1 + k*(sigma0 + sigma1_m * Am[j])));
            // Current forcing
            M3[j] = double(( ((4 * double(5.1419e-10))/(mu*Lx*Ly)) * sin(float_Pi*m[j]*(x_0/Lx)) * sin(float_Pi*n[j]*(y_0/Ly)) ) /
                        (1+k*(sigma0 + sigma1_m * Am[j])));
        }
    }
};


class PERC
{
    std::vector<double> spatial, M1,M2,M3,RC,mem1,mem2,mem3 , MM1, MM2, MM3;
    int MODES,counter,midiNoteNumber,Material;
	float Tension, E,Damping,Thickness, Volume;
    const float pi = 3.141592653589793238462643f, DURATION = 44100;

    ScopedPointer<Membrane>         membrane;
    
public:
	PERC(int m, float vol, float fs, float stiffness,float damp,float thick, int mat,float volume) :
                    midiNoteNumber(m), E(stiffness), Damping(damp),Thickness(thick), Material(mat), Volume(volume)
	{
        float Tension, volDens,  Lx , Ly  , nu = 0.3;
        
        // Map parameters: Material option will define density and E while note number will define the sizes
        //////////////////////////////////////////////
        
        //--- Material
        if (Material == 1)
            {volDens = 3050;}
        else
            {volDens = 10050;}

        //--- Young's Modulus <---> Stiffness! (probably the parameter with the greatest effect...)
        E = (1e11-1e09)*E + 1e09; // map (0.0,1.0) to physical Y.M. values: (1e09,1e11)

        //--- Dimensions
           float midiMax = 108, midiMin = 21, Lx_max = 1.6, Lx_min = 0.22;
            Lx = (Lx_max-Lx_min) * (midiNoteNumber-midiMin) / (midiMax - midiMin) + Lx_min;
            Lx = Lx_max - Lx + Lx_min; Ly = Lx;
            //Ly = Lx + uni(rng);
        
        //--- Constant tension: changing this, does not alter result much
            Tension = 100;
        
        //--- Calculate other parameters and find "allowable" number of modes
        float   D =  (Thickness*Thickness*Thickness*E) / ((12*(1-nu*nu))),
                mu = (Thickness*volDens), sigma0 = Damping, sigma1_m = 0.1*Damping,
                GAMMA =  - 4.0 * 44100.0 * 44100.0, test1 = float_Pi / Lx , test2 = float_Pi / Ly,
                ALPHA = double(D/mu)*double(test1*test1*test1*test1 + (2* float_Pi*float_Pi*float_Pi*float_Pi)/(Lx*Lx*Ly*Ly)+test2*test2*test2*test2),
                BETA =  - double(Tension/mu)*double(pow(pi/Lx,2) + pow(pi/Ly,2)),
                DELTA = BETA  * BETA - 4 * ALPHA * GAMMA;
        float gammaMin = 5, gammaMax = 10;
        float gamma = (gammaMax - gammaMin) * (sigma0) + gammaMin;
        
        int Modes = abs(static_cast<int> (sqrt((sqrt(DELTA) -  BETA) / (2*ALPHA)))) / 2 ; //Bind number of modes so that solution doesn't explode!
        
        // Debugging prints
        //std::cout << "Lx = " << Lx << " , Ly = " << Ly << " , H = " << Thickness << " , D = " << D <<  " , mu = " <<  mu  << std::endl;
        //std::cout << "Modes,E, volDens , Tension " << Modes << "  , "<< E << " ,  " << volDens<<" , " << Tension<< std::endl;
        
        if (Modes > 500) Modes = 500; // Efficiency Fudginess!

        counter = 0;
        MM1.resize(Modes*Modes); MM2.resize(Modes*Modes); MM3.resize(Modes*Modes);
   
    
        membrane = new Membrane(Lx,Ly,Tension,sigma0,sigma1_m,E,Thickness,volDens,Modes,gamma);
	}

	
    // DSP CODE
	void renderToBuffer(float* outbuf, int samples, bool isDown)
	{
        gl_RC = membrane -> RC;
        gl_spatial = membrane->spatial;
        gl_M1 = membrane -> M1;
        gl_M2 = membrane -> M2;
        gl_M3 = membrane -> M3;
        gl_mem2 = MM2;
        gl_mem3 = MM3;
        float gamma = membrane -> gamma;
        int MODES = membrane->Modes; MODES = MODES * MODES;
        
        
        for (int n = 0; n < samples; n++) { //~~~~~~~~~~~ For every sample of the buffer, produce a temporal component vector.
            if (counter > 2) {
                if (counter < 200)  //~~~~~~~~~~~ Here, we will include the forcing samples RC.
                    for (int i=0; i < MODES; i++)
                        gl_mem1[i] = gl_M1[i] * gl_mem2[i] + gl_M2[i]*gl_mem3[i] + gl_M3[i]*gl_RC[counter];
                
                else                //~~~~~~~~~~~ Here, the forcing is over. Only previous states affect the vibration
                    for (int i=0; i < MODES; i++)
                        gl_mem1[i] = gl_M1[i] * gl_mem2[i] + gl_M2[i] * gl_mem3[i];
            }
            

            outbuf[n] = (Volume / 70.0) * double(exp(-gamma*counter/44100))*double( (1/0.3) * 1e15 * std::inner_product(begin(gl_spatial), end(gl_spatial), begin(gl_mem1), 0.0));

            gl_mem3 = gl_mem2; gl_mem2 = gl_mem1; // Update the previous displacement vectors
            
            // After buffer has been "recycled", we need a way to keep an index of the previous states, mem1,mem2,mem3
            MM2 = gl_mem2; MM3 = gl_mem3;
            counter++;
        }
        
    }
    

	bool is_alive() {
        return counter < DURATION;
	}

	int getIndex() {
		return midiNoteNumber;
	}
};

#endif  // PHYSICALMODEL_H_INCLUDED
