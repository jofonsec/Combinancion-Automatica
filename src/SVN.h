#ifndef SingleVertexNeighborhood_h
#define SingleVertexNeighborhood_h

template<class EOT> 
class SingleVertexNeighborhood: public eoMonOp<EOT>{

public:
        /// Constructor principal
        SingleVertexNeighborhood(int _noAnclas, int _geneTotal, int _nodos, double _minimo, double _maximo, eoIncrementorParam<unsigned> * gen):
        noAnclas(_noAnclas), T(_geneTotal), nodos(_nodos), min(_minimo), max(_maximo), nbGenerations(gen){}

        /// Nombre de la clase
        virtual std::string className() const { return "Single vertex neighborhood mutation"; }

        bool operator()(EOT& _Individuo){
                double s = pow(max-min, 2); //Calculo del area de trabajo, es fijo
                double idealDistance = sqrt(s/nodos); // es la distancia ideal, es fijo
                double r; //Es el radio aproximado, se calcula por generacion
                double f, theta=0;

                //Obtiene la generación actual
                int t = atoi(nbGenerations->getValue().c_str());
                r = idealDistance * (1-(t/T));

                for(int i=noAnclas; i<nodos; i++){
                        theta = eo::rng.uniform(0.0,360.0);
                        _Individuo[i*2] = _Individuo[i*2] + r * cos(theta);
                        _Individuo[i*2+1] = _Individuo[i*2+1] + r * sin(theta);
                }
                return true;
        }

private:
        double min, max;
        int T; //total de generaciones, es fijo
        int noAnclas;
        int nodos; //Total de nodos, es fijo
        eoIncrementorParam<unsigned>  * nbGenerations; /*!< numero de generaciones */

};

#endif