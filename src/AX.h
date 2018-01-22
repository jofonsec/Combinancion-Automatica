#ifndef ArithmeticCrossover_H_
#define ArithmeticCrossover_H_

#include <algorithm>    // swap_ranges
#include <utils/eoRNG.h>
#include <es/eoReal.h>
#include <utils/eoRealVectorBounds.h>

template <class EOT>
class ArithmeticCrossover : public eoQuadOp<EOT>
{
public:
   ArithmeticCrossover(const double& _alpha = 0.0) :
   bounds(eoDummyVectorNoBounds), alpha(_alpha), range(1+2*_alpha) {}
 
   ArithmeticCrossover(eoRealVectorBounds & _bounds, const double& _alpha = 0.0) :
   bounds(_bounds), alpha(_alpha), range(1+2*_alpha) {}
 
   virtual std::string className() const { return "ArithmeticCrossover"; }
 
bool operator()(EOT& _eo1, EOT& _eo2){
unsigned i;
double r1, r2, fact;
double alphaMin = -alpha;
double alphaMax = 1+alpha;
   if (alpha == 0.0)            // no check to perform
      fact = -alpha + rng.uniform(range); // in [-alpha,1+alpha)
   else {                      // look for the bounds for fact
        for (i=0; i<_eo1.size(); i++){
                r1=_eo1[i];
                r2=_eo2[i];
                if (r1 != r2) {      // otherwise you'll get NAN's
                        double rmin = std::min(r1, r2);
                        double rmax = std::max(r1, r2);
                        double length = rmax - rmin;
                        if (bounds.isMinBounded(i)){
                                alphaMin = std::max(alphaMin, (bounds.minimum(i)-rmin)/length);
                                alphaMax = std::min(alphaMax, (rmax-bounds.minimum(i))/length);
                        }
                        if (bounds.isMaxBounded(i)){
                                alphaMax = std::min(alphaMax, (bounds.maximum(i)-rmin)/length);
                                alphaMin = std::max(alphaMin, (rmax-bounds.maximum(i))/length);
                        }
                }
        }
        fact = alphaMin + (alphaMax-alphaMin)*rng.uniform();
        }

        for (i=0; i<_eo1.size(); i++){
                r1=_eo1[i];
                r2=_eo2[i];
                _eo1[i] = fact * r1 + (1-fact) * r2;
                _eo2[i] = (1-fact) * r1 + fact * r2;
                }
        return true;         // shoudl test if fact was 0 or 1 :-)))
}

protected:
        eoRealVectorBounds & bounds;
        double alpha;
        double range;                    // == 1+2*alpha
};

#endif // ArithmeticCrossover_H
