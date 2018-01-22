#include <Individuo.h>
#include <eoEvalFunc.h>

class localizacionEvalPenal : public eoEvalFunc<Individuo>
{
    double DisReal[500][500];
    double vecAnclas[500];
    double radio;
    int NoAnclas; //Cantidad nodos anclas
    int nodos; //Total de nodos

public:
    void operator()(Individuo & _Individuo);
    void guardarDisReal(double _DisReal[500][500], int _NoAnclas, int _nodos, double _radio);
    void guardarAnclas(double _vecAnclas[500]);
    void muestraDisReal();
    void muestraAnclas();
};
