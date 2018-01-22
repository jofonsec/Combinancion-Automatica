#ifndef SwapMutation_h
#define SwapMutation_h


template<class Chrom> 
class SwapMutation: public eoMonOp<Chrom>
{
public:
        SwapMutation(const unsigned _howManySwaps=1): howManySwaps(_howManySwaps){
                // consistency check
                if(howManySwaps < 1)
                        throw std::runtime_error("Invalid number of swaps in eoSwapMutation");
        }

        virtual std::string className() const { return "eoSwapMutation"; }

        bool operator()(Chrom& chrom){
                unsigned i, j;
                
                for(unsigned int swap = 0; swap < howManySwaps; swap++){
                        // generate two different indices
                        i=eo::rng.random(chrom.size());
                        do j = eo::rng.random(chrom.size()); while (i == j);
                        // swap
                        std::swap(chrom[i],chrom[j]);
                }
                return true;
        }

private:
        unsigned int howManySwaps;

};

#endif