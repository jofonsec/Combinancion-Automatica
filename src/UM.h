#ifndef UniformMutation_h
#define UniformMutation_h

#include <algorithm>    // swap_ranges
#include <utils/eoRNG.h>
#include <es/eoReal.h>
#include <utils/eoRealVectorBounds.h>

template<class EOT> 
class UniformMutation: public eoMonOp<EOT>{

public:
        UniformMutation(const double& _epsilon, const double& _p_change = 1.0): 
        homogeneous(true), bounds(eoDummyVectorNoBounds), epsilon(1, _epsilon),
        p_change(1, _p_change) {}

        UniformMutation(eoRealVectorBounds & _bounds, const double& _epsilon, const double& _p_change = 1.0):
        homogeneous(false), bounds(_bounds), epsilon(_bounds.size(), _epsilon),
        p_change(_bounds.size(), _p_change){
        // scale to the range - if any
        for (unsigned i=0; i<bounds.size(); i++)
                if (bounds.isBounded(i))
                        epsilon[i] *= _epsilon*bounds.range(i);
        }

        UniformMutation(eoRealVectorBounds & _bounds, const std::vector<double>& _epsilon, const std::vector<double>& _p_change):
        homogeneous(false), bounds(_bounds), epsilon(_epsilon),
        p_change(_p_change) {}

        virtual std::string className() const { return "UniformMutation"; }

        bool operator()(EOT& _eo){
                bool hasChanged=false;
                if (homogeneous)             // No hay objeto de limites "bounds"
                for (unsigned i=0; i<_eo.size(); i++){
                        if (rng.flip(p_change[0])){
                                _eo[i] += 2*epsilon[0]*rng.uniform()-epsilon[0];
                                hasChanged = true;
                        }
                }
                else{
                        // sanity check ?
                        if (_eo.size() != bounds.size())
                                throw std::runtime_error("Invalid size of indi in UniformMutation");

                        for (unsigned i=0; i<_eo.size(); i++)
                                if (rng.flip(p_change[i])){
                                        // check the bounds
                                        double emin = _eo[i]-epsilon[i];
                                        double emax = _eo[i]+epsilon[i];
                                        if (bounds.isMinBounded(i))
                                                emin = std::max(bounds.minimum(i), emin);
                                        if (bounds.isMaxBounded(i))
                                                emax = std::min(bounds.maximum(i), emax);
                                        _eo[i] = emin + (emax-emin)*rng.uniform();
                                        hasChanged = true;
                                }
                }
                return hasChanged;
        }
 
private:
        bool homogeneous;   // == no bounds passed in the ctor
        eoRealVectorBounds & bounds;
        std::vector<double> epsilon;     // the ranges for mutation
        std::vector<double> p_change;    // the proba that each variable is modified
};

#endif