#include <sstream>

#include<vector>
#include<iterator>
#include<string>
#include<iostream>

#include "GnuplotMonitor.h"

#include "gnuplot-iostream.h"

#include <unistd.h> //API POSIX

eoMonitor& GnuplotMonitor::operator() (void)
{
    //Actualiza el archivo eoFileMonitor
    //eoFileMonitor::operator()();

    //Obtener anclas desde archivo txt
    std::ifstream fin,fin1,fin2,fin3,fin4;

    int generacion = atoi(nbGenerations->getValue().c_str());

	std::vector<double> fitness(generacion), fitness1(generacion), fitness2(generacion), fitness3(generacion);
	std::vector<double> media(generacion),media1(generacion),media2(generacion),media3(generacion);
	std::vector<double> Diversity(generacion), Diversity1(generacion),Diversity2(generacion),Diversity3(generacion);
    //std::cout << "Press Ctrl-C to quit (closing gnuplot window doesn't quit)." << std::endl;

    if(activate==1){
        //Obtener anclas desde archivo txt
        fin.open(filename);

        if(islas==1){
            fin1.open("Resumen1.xg"); fin2.open("Resumen2.xg"); fin3.open("Resumen3.xg"); fin4.open("Resumen4.xg");

            if(!fin1.is_open() || !fin2.is_open()|| !fin3.is_open()||!fin4.is_open()){
            std::cout << "Error mientras se intentaba leer el archivo "<< "Resumen1.xg, 2, 3 y 4" << ". Por favor verificar que existan!" << std::endl;
            exit(1);
            }
            else {
                for (int i =0 ; i < generacion ; i++){
                    fin1 >> fitness[i];
                    fin2 >> fitness1[i] >> media1[i] >> Diversity1[i];
                    fin3 >> fitness2[i] >> media2[i] >> Diversity2[i];
                    fin4 >> fitness3[i];
                }
            fin1.close(); fin2.close(); fin3.close(); fin4.close();
            }
        gp <<"set autoscale\n";
        gp << "plot '-' with lines title 'Isla1', '-' with lines title 'Isla2', '-' with lines title 'Isla3', '-' with lines title 'Isla4'\n";
        gp.send1d(fitness);
        gp.send1d(fitness1);
        gp.send1d(fitness2);
        gp.send1d(fitness3);
        }
        else{
            if(!fin.is_open()){
            std::cout << "Error mientras se intentaba leer el archivo "<< filename << ". Por favor verificar que exista!" << std::endl;
            exit(1);
            }
            else {
                for (int i =0 ; i < generacion ; i++)
                fin >> fitness[i] >> media[i] >> Diversity[i];
            fin.close();
            }
        gp <<"set autoscale\n";
        gp << "plot '-' with lines title 'Fitness', '-' with lines title 'Media'\n";
        gp.send1d(fitness);
        gp.send1d(media);
        }
    }
}

void GnuplotMonitor::setGen(eoIncrementorParam<unsigned> * gen){
    nbGenerations = gen;
}
