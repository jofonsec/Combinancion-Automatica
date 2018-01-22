//Se debe cambiar el nombre de uno de los archivos .sav creado por el siguiente:
// ContarLineas.sav

#include "graphError.h"

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <C:/Program\ Files/boost/boost_1_62_0/boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"

graphError::graphError(std::string _filename, int setGeneracion, int numberGeneration, int _nodos, int _noAnclas, double _max)
{
    nodos = _nodos;
    filename = _filename;
    noAnclas = _noAnclas;
    int margen = 3;

    if (_max>25) margen=6;

    //Cargar posiciones
    cargar_unknow();

    //Ejecucion de comandos en el shell
    std::string archivo2=filename+std::to_string(setGeneracion)+ext;
    std::string s="cp "+archivo2+" ContarLineas.sav";
    char *comando = strdup(s.c_str());
    system(comando);
    system("sed -i '1,2d' *.sav");

    for(int k=setGeneracion;k<numberGeneration+1;k+=setGeneracion){

        archivo=filename+std::to_string(k)+ext; //Define el nombre del archivo que guarda la poblacion
        filas = contador_filas(); //Se debe cambiar el nombre del archivo desde el cual se contaran las lineas
        cargar_poblacion(archivo);
        ordenar_matriz();
        mejor_individuo();

        //Graficar
        Gnuplot gp;
        std::vector<std::vector<std::pair<double, double> > > all_segments;
        std::vector<std::pair<double, double> > xy_pts_anclas;
        std::vector<std::pair<double, double> > xy_pts_unknow;
        //Grafica los segmentos (Error)
        int contador=0; //Contador para segmentos
        for(int j=0; j<nodos-noAnclas; j++) {
            std::vector<std::pair<double, double> > segment;
            for(int i=0; i<1; i++) {
                segment.push_back(std::make_pair(unknowX[contador], unknowY[contador]));
                segment.push_back(std::make_pair(individuoX[contador], individuoY[contador]));
                contador++;
            }
            all_segments.push_back(segment);
        }

        //Grafica los puntos
        for(int i=0; i<noAnclas; i++) xy_pts_anclas.push_back(std::make_pair(anclasX[i], anclasY[i]));

        for(int i=0; i<nodos-noAnclas; i++) xy_pts_unknow.push_back(std::make_pair(unknowX[i], unknowY[i]));

        gp << "set yrange [0:"<<_max+margen<<"]\n";
        //gp << "set autoscale\n";
        gp << "set title 'Generaci�n "<<std::to_string(k)<<"' font 'Bookman,20'\n";
        gp << "set encoding iso_8859_1\n"; //Define los acentos para postscript
        gp << "set grid\n";
        gp << "set xlabel 'metros'\n";
        gp << "set ylabel 'metros'\n";
        gp << "set key font 'Helvetica,18'\n";//Modifica la leyenda
        gp << "set terminal postscript enhanced color font 'Helvetica,15'\n";//Cambia de terminal
        gp << "set output '| ps2pdf - Generacion_"<<std::to_string(k)<<".pdf'\n";//guarda en pdf
        gp << "set style line 1 lc rgb '#00DD00' pt 7 ps 1.5 \n";
        gp << "set style line 2 lc rgb 'blue' pt 6 ps 1.5 \n";
        gp << "plot '-' with line title 'Error', '-' with points ls 1 title 'Nodos Anclas',"
        << "'-' with points ls 2 title 'Posici�n Real'\n";
        //<< "'-' with vectors title 'pts_B'\n";
        gp.send2d(all_segments);
        gp.send1d(xy_pts_anclas);
        gp.send1d(xy_pts_unknow);
    }
}

void graphError::mejor_individuo(){
    for(int i=0;i<nodos-noAnclas;i++){
         individuoX[i]=Matriz[0][i*2+(noAnclas*2+2)];
         individuoY[i]=Matriz[0][i*2+(noAnclas*2+3)];
    }
}

void graphError::ordenar_matriz(){
    double temp;
    for(int k=0;k<filas-2;k++){
        for(int i=0;i<filas-2;i++){
            if(Matriz[i][0]>Matriz[i+1][0]){
                for(int j=0;j<nodos*2+2;j++){
                    temp = Matriz[i][j];
                    Matriz[i][j] = Matriz[i+1][j];
                    Matriz[i+1][j]=temp;
                }
            }
        }
    }
}

void graphError::cargar_unknow(){
    std::fstream fin;
    fin.open(archivo1.c_str());

    //Mensaje de error si no encuentra el archivo Allposition.txt
    if(!fin.is_open()){
        std::cout << "Error mientras se intentaba leer el archivo "<< "Allposition.txt" << ". Por favor verificar que exista!" << std::endl;
        exit(1);
    }
    else {
        for(int j=0;j<noAnclas;j++){
            fin >> anclasX[j] >> anclasY[j];
        }

        for(int j=0;j<nodos-noAnclas;j++){
            fin >> unknowX[j] >> unknowY[j];
        }
        fin.close();
    }
}

void graphError::cargar_poblacion(std::string archivo){
    std::fstream fin;
    fin.open(archivo.c_str());
    for(int i=0;i<filas;i++)
        for(int j=0;j<nodos*2+2;j++){
            fin >> Matriz[i][j];
        }
    fin.close();
}

int graphError::contador_filas(){
    int contador=0, contar=0, hc;
    //Contar la cantidad de lineas y posicion de un archivo, stdio
    FILE *leer = fopen("ContarLineas.sav", "r");
    while (EOF != (hc=getc(leer))){
        contador++;
        if ('\n'==hc){
            contar++;
        }
    }
    return contar;
}
