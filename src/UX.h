#ifndef UniformCrossover_H_
#define UniformCrossover_H_

#include <algorithm>    // swap_ranges
#include <utils/eoRNG.h>
#include <es/eoReal.h>
#include <utils/eoRealVectorBounds.h>

template<class EOT> 
class UniformCrossover: public eoQuadOp<EOT>
{
        public:
                UniformCrossover(const float& _preference = 0.5): preference(_preference){
                        if ( (_preference <= 0.0) || (_preference >= 1.0) )
                                std::runtime_error("UniformCrossover --> invalid preference");
                }

                virtual std::string className() const { return "UniformCrossover"; }

                bool operator()(EOT& _eo1, EOT& _eo2){
                        if ( _eo1.size() != _eo2.size())
                                std::runtime_error("UniformCrossover --> Tamaño chromosomas no coinciden" );
                        bool changed = false;
                        for (unsigned int i=0; i<_eo1.size(); i++){
                                if (rng.flip(preference))//lanzamiento de moneda sesgada (probabilidad de cara o cruz son diferentes). Si preference es mayor hay altas probabilidades de que se cambien todos los genes
                                        if (_eo1[i] != _eo2[i]){
                                                double tmp = _eo1[i];
                                                _eo1[i]=_eo2[i];
                                                _eo2[i] = tmp;
                                                changed = true;
                                        }
                        }
                        return changed;
                }
        private:
                float preference;
};

#endif // UniformCrossover_H
