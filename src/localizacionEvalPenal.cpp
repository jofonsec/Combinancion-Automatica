#include <localizacionEvalPenal.h>
using namespace std;

void localizacionEvalPenal::operator()(Individuo & _Individuo){

double ErrorTotal = 0.0, DisEstimada = 0.0;
double Error1 = 0.0, Error2 = 0.0;
int Error3 = 0;
int contador =0;

    for (int i= 0; i< NoAnclas; i++){//Primero se evalua con la posicion de los nodos anclas
        for (int j= NoAnclas; j < _Individuo.size()/2; j++){
            if (DisReal[i][j] != 0){//Primero se evalua los nodos que se sabe estan conectados
                DisEstimada = sqrt(pow(vecAnclas[i*2] - _Individuo[j*2],2) + pow(vecAnclas[i*2+1] - _Individuo[j*2+1],2));
                Error1 = abs (DisEstimada - DisReal[i][j]) + Error1;//Error acumulado con respecto a los nodos anclas
            }
            else{
                DisEstimada = sqrt(pow(vecAnclas[i*2] - _Individuo[j*2],2) + pow(vecAnclas[i*2+1] - _Individuo[j*2+1],2));
                if(DisEstimada <= radio)Error3++;//Error acumulado de falsas conexiones
            }
        }
    }

    for (int i= NoAnclas; i<_Individuo.size()/2;i++){
        contador++;//Permite ir desplazando la columna
        for(int j= NoAnclas+contador; j< _Individuo.size()/2; j++){
            if (DisReal[i][j] != 0){//Comparacion entre nodos con posiciones desconocidas
                DisEstimada = sqrt(pow(_Individuo[i*2] - _Individuo[j*2],2) + pow(_Individuo[i*2+1] - _Individuo[j*2+1],2));
                Error2 = abs(DisEstimada - DisReal[i][j]) + Error2;//Error acumulado de nodos desconocidos
            }
            else{
                DisEstimada = sqrt(pow(_Individuo[i*2] - _Individuo[j*2],2) + pow(_Individuo[i*2+1] - _Individuo[j*2+1],2));
                if(DisEstimada <= radio) Error3++;//Error acumulado de falsas conexiones
            }
        }
    }
    // Error2 0<==>1200 and Error1 0<==>200 and Error3 0<==>600
    ErrorTotal = Error2/10 + 3 * Error1 + 8 * Error3;

    // Finalmente se debe agregar el fitness al cromosoma
    _Individuo.fitness(ErrorTotal);
}

void localizacionEvalPenal::guardarDisReal(double _DisReal[500][500], int _NoAnclas, int _nodos, double _radio){
    nodos = _nodos;
    NoAnclas = _NoAnclas;
    radio = _radio;

    for (int i=0; i<nodos ; i++){
        for (int j=0; j<nodos; j++){
            DisReal[i][j] = _DisReal[i][j];
        }
    }
}

void localizacionEvalPenal::guardarAnclas(double _vecAnclas[500]){

     for (int i=0 ; i<NoAnclas*2 ; i++){
        vecAnclas[i] = _vecAnclas[i];
     }

}

void localizacionEvalPenal::muestraDisReal(){
   for (int i=0; i < nodos; i++){
        for (int j=0; j < nodos ; j++){
            std::cout << DisReal[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}

void localizacionEvalPenal::muestraAnclas(){
       for (int i=0; i < NoAnclas*2; i++){
            std::cout << vecAnclas[i] << "  ";
        }
}
