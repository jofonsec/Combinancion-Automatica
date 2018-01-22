#ifndef SimulatedBinaryCrossover_H_
#define SimulatedBinaryCrossover_H_

#include <algorithm> 
#include <utils/eoParser.h>
#include <utils/eoRNG.h>
#include <es/eoReal.h>
#include <utils/eoRealBounds.h>
#include <utils/eoRealVectorBounds.h>

template <class EOT>
class SimulatedBinaryCrossover : public eoQuadOp<EOT>
{
public:
        SimulatedBinaryCrossover(const double& _eta = 1.0) :
        bounds(eoDummyVectorNoBounds), eta(_eta), range(1) {}

        SimulatedBinaryCrossover(eoRealVectorBounds & _bounds, const double& _eta = 1.0, unsigned int& _NoAnclas = 10) :
        bounds(_bounds), eta(_eta), range(1), NoAnclas(_NoAnclas) {}

        SimulatedBinaryCrossover(eoParser & _parser) :
        bounds (_parser.getORcreateParam(eoDummyVectorNoBounds, "objectBounds", "Bounds for variables", 'B', "Variation Operators").value()) , eta (_parser.getORcreateParam(1.0, "eta", "SBX eta parameter", '\0', "Variation Operators").value()) ,
        range(1) {}

        bool operator()(EOT& _individuo1, EOT& _individuo2){
                double r1, r2, beta;
                
                for (unsigned i=NoAnclas*2; i<_individuo1.size(); i++){
                        
                        double u = rng.uniform(range);
                        if ( u <= 0.5 )
                                beta = exp( (1/(eta+1))*log(2*u));
                        else
                                beta = exp((1/(eta+1))*log(1/(2*(1-u))));

                        r1=_individuo1[i];
                        r2=_individuo2[i];
        
                        _individuo1[i] =0.5*((1+beta)*r1+(1-beta)*r2);
                        _individuo2[i] =0.5*((1-beta)*r1+(1+beta)*r2);

                        if(!(bounds.isInBounds(i,_individuo1[i]))) //Verifica si elemento se encuentra dentro de los limites especificados
                                bounds.foldsInBounds(i,_individuo1[i]); //Ingresa un valor dentro de los limites
                        if(!(bounds.isInBounds(i,_individuo2[i])))
                                bounds.foldsInBounds(i,_individuo2[i]);
                }
        }

protected:
        eoRealVectorBounds & bounds;
        double eta;
        double range;
        int NoAnclas;
        int contador;                    // == 1
};

#endif // SimulatedBinaryCrossover_H
