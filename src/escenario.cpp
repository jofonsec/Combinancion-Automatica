#include "escenario.h"

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <random>

escenario::escenario(double _minimo, double _maximo, double _radio, int _nodos, int _aleatorio, int _NoAnclas)
{
    //ctor
    std::mt19937 generador;
    std::normal_distribution<double> normal(0,1);
    nodos=_nodos;
    NoAnclas=_NoAnclas;
    double f, f1;
    double distancia=0.0;
    srand(_aleatorio);

    //Generar las coordenadas de forma aleatoria
    for (int i=0; i<_nodos; i++){
        f = (double)rand() / RAND_MAX;
        f1 = (double)rand() / RAND_MAX;
        xt.push_back(_minimo + f * (_maximo-_minimo));
        yt.push_back(_minimo + f1 * (_maximo-_minimo));
    }

    //Generar vector de coordenadas nodos anclas
    for (unsigned i=0 ; i < _NoAnclas ; i++){
        anclas[i*2] = xt[i];
        anclas[i*2+1] = yt[i];
    }

    //Matrix of zeros
    for (int i=0; i<_nodos; i++){
        for (int j=0; j<_nodos ; j++){
        DisRSSI[i][j]=0;
        DisReal[i][j]=0;
        }
    }

    //Llenar con distancia real
    for (int i=0; i<_nodos; i++){
        for (int j=_NoAnclas; j<_nodos ; j++){
            if(j>i){
                distancia = sqrt(pow(xt[i] - xt[j],2) + pow(yt[i] - yt[j],2));
                //std::cout <<" i " << i <<" j " << j <<" : " << xt[i] <<" y " << yt[i] << std::endl;
                if (distancia < _radio){
                    //Aplicar ruido gaussiano a las distancias
                    DisRSSI[i][j]= distancia*(1.0+0.1*normal(generador));
                    DisReal[i][j]= distancia;
                }
            }
        }
    }
}

escenario::escenario(int _nodos, int _NoAnclas){
    nodos=_nodos;
    NoAnclas=_NoAnclas;
    //Obtener anclas desde archivo txt
    std::ifstream fin;
    double data1,data2,data3;
    fin.open("Anclas.txt");

    if(!fin.is_open()){
        std::cout << "Error mientras se intentaba leer el archivo "<< "Anclas.txt" << ". Por favor verificar que exista!" << std::endl;
        exit(1);
    }
    else {
        for (int i =0 ; i < _NoAnclas; i++)
            fin >> anclas[i*2] >> anclas[i*2+1];
    fin.close();
    }

    //Obtener matriz de distancia desde txt
    fin.open("Matriz.txt");

    if(!fin.is_open()){
        std::cout << "Error mientras se intentaba leer el archivo "<< "Matriz.txt" << ". Por favor verificar que exista!" << std::endl;
        exit(1);
    }
    else {
        for (int i = 0 ; i < _nodos; i++){
            for(int j = 0; j < _nodos; j++)
            fin >> DisRSSI[i][j];
        }
    fin.close();
    }
}

void escenario::guardarNodos(){
    std::fstream dataarch;
    dataarch.open("Allposition.txt", std::ios::out);

    for(int i=0;i<nodos;i++){
    dataarch<< xt[i] <<"  "<< yt[i]<< std::endl;
    }
    dataarch.close();
}

void escenario::guardarAnclas(){
    std::fstream dataarch;
    dataarch.open("Anclasposition.txt", std::ios::out);

    for(int i=0;i<NoAnclas;i++){
    dataarch<< anclas[i*2] <<"  "<< anclas[i*2+1]<< std::endl;
    }
    dataarch.close();
}

void escenario::guardarDisRSSI(){
    std::fstream dataMa;
    dataMa.open("MatrizRSSI.txt", std::ios::out);

    for (int i=0; i < nodos; i++){
        for (int j=0; j < nodos ; j++){
            dataMa << DisRSSI[i][j] << "  ";
        }
        dataMa << std::endl;
    }
    dataMa.close();
}

void escenario::guardarDisReal(){
    std::fstream dataMa1;
    dataMa1.open("MatrizReal.txt", std::ios::out);

    for (int i=0; i < nodos; i++){
        for (int j=0; j < nodos ; j++){
            dataMa1 << DisReal[i][j] << "  ";
        }
        dataMa1 << std::endl;
    }
    dataMa1.close();
}

double escenario::obtenerDisRSSI(int i, int j){
    double valor = DisRSSI[i][j];
    return valor;
}

double escenario::obtenerAnclas(int i){
    double valor = anclas[i];
    return valor;
}

double escenario::obtenerPX(int i){
    double valorX = xt[i];
    return valorX;
}
double escenario::obtenerPY(int i){
    double valorY = yt[i];
    return valorY;
}
escenario::escenario(double _radio, int _nodos, int _NoAnclas){
    nodos=_nodos;
    double distancia=0.0;

    //Cargar escenario de prueba
    //Iniciar los vectores de prueba
    double nodosX[]={0.184942, 16.1406, 15.3417, 18.9813, 2.68624, 16.2609, 15.2098, 9.68169, 11.3004, 12.3374,17.1685,15.5278,13.2145,
    11.0898,19.9731,9.85504,9.24223,6.86728,16.0448,13.9824,17.0251,10.4764,10.922,12.2861,6.73605,10.5081,15.2141,19.7803,9.28373,7.1688,
    9.05362,2.19367,17.0147,18.6059,19.0521,8.55068,16.9158,19.7375,7.74804,16.245,16.3182,13.5789,15.5669,8.62026,1.10172,15.0731,15.2672,
    5.04959,6.39668,7.41478,18.8995,3.51024,19.7845,10.8853,17.2564,15.6652,17.0061,4.17249,11.8247,14.4371,10.1596,18.2061,19.2425,14.0916,
    8.40968,2.47932,8.24854,2.81198,11.1844,3.88501,9.13724,7.77734,13.881,18.0358,15.7115,12.7427,6.07074,13.0955,19.7351,8.11426,2.19916,
    6.1092,17.8857,13.5313,4.7908,11.8284,11.2101,12.9398,0.151982,4.04004,17.0763,2.17475,6.57369,9.39604,13.6235,9.03653,15.0359,13.7382,
    16.7211,5.86688};
    xt.insert(xt.begin(), nodosX, nodosX+100);

    double nodosY[]={16.0967 ,11.1411 ,1.33183 ,14.9858 ,14.7337 ,19.3292 ,18.197 ,6.8862 ,6.95822 ,11.984 ,4.81887 ,2.60323 ,7.53685 ,
    15.9008 ,7.50755 ,5.32853 ,2.85775 ,7.48497 ,4.0791 ,6.81234 ,1.3245 ,6.27216 ,5.49638 ,2.12653 ,6.24531 ,10.9391 ,0.00976592 ,6.93258 ,
    11.2223 ,17.821 ,19.7986 ,17.2448 ,18.0462 ,1.92816 ,10.1096 ,11.6971 ,17.0812 ,0.46144 ,0.899075 ,16.3103 ,10.1621 ,18.8263 ,0.976592 ,
    11.7539 ,15.9795 ,19.1485 ,14.6385 ,3.11533 ,17.973 ,11.2918 ,0.813013 ,12.1647 ,14.7215 ,6.3741 ,14.7038 ,12.1177 ,12.8343 ,12.7781 ,
    11.9224 ,4.90432 ,2.33894 ,14.6715 ,15.1372 ,0.935087 ,14.5781 ,13.0222 ,3.20811 ,1.36357 ,2.59041 ,9.11222 ,2.8309 ,2.85409 ,4.01257 ,
    10.5539 ,2.67647 ,4.84634 ,6.46016 ,19.8071 ,10.0302 ,19.3225 ,15.4241 ,5.19059 ,9.935 ,3.06955 ,4.65896 ,14.4041 ,6.83004 ,12.4516 ,
    9.07376 ,10.5387 ,8.69289 ,0.665914 ,0.123905 ,10.5936 ,13.6784 ,17.5188 ,3.70006 ,6.49312 ,19.2615 ,11.6819};
    yt.insert(yt.begin(), nodosY, nodosY+100);

    //Generar vector de coordenadas nodos anclas
    for (unsigned i=0 ; i < _NoAnclas ; i++){
        anclas[i*2] = xt[i];
        anclas[i*2+1] = yt[i];
    }

    //Matrix of zeros
    for (int i=0; i<_nodos; i++){
        for (int j=0; j<_nodos ; j++){
        DisReal[i][j]=0;
        }
    }

    //Llenar con distancia real
    for (int i=0; i<_nodos; i++){
        for (int j=_NoAnclas; j<_nodos ; j++){
            if(j>i){
                distancia = sqrt(pow(xt[i] - xt[j],2) + pow(yt[i] - yt[j],2));
                if (distancia <= _radio)
                    DisReal[i][j]= distancia;
            }
        }
    }
}
