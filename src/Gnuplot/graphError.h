#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <C:/Program Files/boost/boost_1_62_0/boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"
class graphError
{
    private:
        double Matriz[500][500];
        int filas;
        int nodos;
        int noAnclas;
        std::string filename;
        std::string ext=".sav", archivo , archivo1="Allposition.txt";
        double anclasX[100],anclasY[100], unknowX[1000],unknowY[1000];
        double individuoX[2000], individuoY[2000]; //Almacena el mejor individuo
        //Funciones privadas
        int contador_filas();
        void cargar_poblacion(std::string archivo);
        void cargar_unknow();
        void ordenar_matriz();
        void mejor_individuo();
    public:
        //Grafica el error y posiciones
        graphError(std::string _filename, int setGeneracion, int numberGeneration, int _nodos, int _noAnclas, double _max); //constructor de la clase

};

