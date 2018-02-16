#ifndef evalGAs_h
#define evalGAs_h

#include <eo>
#include <Individuo.h>
#include <IndiInit.h>
#include <localizacionEvalPenal.h>
#include <sys/time.h>

//Operadores de cruce
#include <UX.h>
#include <AX.h>
#include <SBX.h>

//Operadores de mutacion
#include <SWAP.h>
#include <UM.h>
#include <SVN.h>

typedef SimulatedBinaryCrossover<Individuo> SBX;
typedef ArithmeticCrossover<Individuo> AX;
typedef UniformCrossover<Individuo> UX;
typedef SwapMutation<Individuo> SWAP;
typedef SingleVertexNeighborhood<Individuo> SVN;
typedef UniformMutation<Individuo> UM;
typedef eoIncrementorParam<unsigned> Incrementor;

template <class EOT>
class evalGAs : public eoEvalFunc<EOT>{

  protected:
  //Datos necesarios del escenario de prueba
    double _min = 0.0;
    double _max = 200.0;
    unsigned int NoAnclas = 20;
    unsigned int nodos = 120;
    double radio = 40.0;

    double DisReal[500][500];
    double vecAnclas[40];

  public :
    evalGAs (double _DisReal[500][500], double _vecAnclas[500]){
      for (int i=0; i<nodos ; i++){
          for (int j=0; j<nodos; j++){
              DisReal[i][j] = _DisReal[i][j];
          }
      }

      for (int i=0 ; i<NoAnclas*2 ; i++){
         vecAnclas[i] = _vecAnclas[i];
      }
    }
    ~evalGAs(){
      std::cout << "Se ejecuta el destructor" << std::endl;
    }
    void operator()(EOT& _eo){

      /**Conversion de Genotipo a valores de probabilidad**/
        double probabilidades[8];
        double Divisor1 = 0.0, Divisor2 = 0.0, Divisor3 = 0.0;

        for (unsigned int i= 0; i< _eo.size(); i+=7){
          probabilidades[i/7]=_eo[i]*64+_eo[i+1]*32+_eo[i+2]*16+_eo[i+3]*8+_eo[i+4]*4+_eo[i+5]*2+_eo[i+6]*1;
        }

        Divisor1 = probabilidades[0] + probabilidades[1] + probabilidades[2];
        Divisor2 = probabilidades[3] + probabilidades[4] + probabilidades[5];
        Divisor3 = probabilidades[6] + probabilidades[7];

        double PcruzaA = probabilidades[0]/Divisor1;
        double PcruzaB = probabilidades[1]/Divisor1;
        double PcruzaC = probabilidades[2]/Divisor1;
        double PmutationA = probabilidades[3]/Divisor2;
        double PmutationB = probabilidades[4]/Divisor2;
        double PmutationC = probabilidades[5]/Divisor2;
        double PcruzaT = probabilidades[6]/Divisor3;
        double PmutationT = probabilidades[7]/Divisor3;

      /**Parametros de Parser**/
      //Primero se debe definir un parser que lee desde la linea de comandos o un archivo
          //eoParser parser(argc, argv);

      //Se definen los parametros, se leen desde el parser y le asigna el valor
        unsigned seed = time(0);

      //Configuracion parametros algoritmo
          unsigned int POP_SIZE = 50;
          unsigned int numberGeneration = 1500;
          unsigned int Nc = 2;
          double Alpha = 0.5;
          float preferencia = 0.5;
          double epsilon = 5;
          double P_change = 0.79;
//          double PmutationA = 0.33 ; //parser.createParam((double)(0.85), "Pmutacion", "Probabilidad de mutacion de la encapsulacion de SVN y Swap",'Y',"Parametros Algoritmo").value();
//          double PmutationB = 0.33 ; //parser.createParam((double)(0.85), "Pmutacion1", "Probabilidad de mutacion de SVN",'Z',"Parametros Algoritmo").value();
//          double PmutationC = 0.33 ; //parser.createParam((double)(0.5), "Pmutacion2", "Probabilidad de mutacion de Swap",'W',"Parametros Algoritmo").value();
//          double PmutationT = 0.33 ;
//          double PcruzaA = 0.33 ;
//          double PcruzaB = 0.33 ;
//          double PcruzaC = 0.33 ;
//          double PcruzaT = 0.33 ;
          double sizeTorneo = 8;
          double sizeElist = 2;
          double sizeTorneo1 = 2;
      //Parametro de tiempo
          struct timeval ti, tf;
          double tiempo;

std::cout<<"Comienza la inicializacion"<< std::endl;

      /**COMIENZO DEL ALGORITMO**/
          //Define la representaciOn (Individuo)
              Individuo cromosoma;

          //Generar una subclase de la clase de la funcion de evaluacion
              localizacionEvalPenal Fitness;

          //Criterio de parada
              eoGenContinue<Individuo> parada(numberGeneration);

          /** CRUZA **/

              // Generar los limites para cada gen
              std::vector<double> min_b;
              std::vector<double> max_b;
              for(int i=0; i<nodos*2; i++) {
                      min_b.push_back(_min);
                      max_b.push_back(_max);
                  }
              eoRealVectorBounds bounds(min_b, max_b);

          /** CRUZAMIENTO**/
              //Inicializar operador de cruce SBX
              SBX crossoverA(bounds, Nc, NoAnclas);

              //Inicializa operador de cruce aritmetico
              AX crossoverB(bounds,Alpha);

              //Inicializa operador de cruce uniforme
              UX crossoverC(preferencia);

              eoPropCombinedQuadOp<Individuo> crossover(crossoverA,PcruzaA);
              crossover.add(crossoverB,PcruzaB);
              crossover.add(crossoverC,PcruzaC);

          /** MUTACION **/
              //Subclase de mutacion paper IEEE
                  //Se inicializa el contador de generaciones
                  //Incrementor generationCounter("Gen.");
                  Incrementor *generationCounter = new Incrementor("Gen.");

              //Operador especifico para el problema
              SVN mutationA(NoAnclas, numberGeneration, nodos, _min, _max, & *generationCounter);
                                                                          //el uso del ampersand permite compartir la direccion en memoria de una variable.
              //Mutacion incluida en EO, permite llegar mas rapido a un fitness de 600
              SWAP mutationB;

              //mutationA.setGen(& generationCounter);
              UM mutationC(bounds, epsilon, P_change);

              //Combina operadores de mutacion con su respectivo peso
              eoPropCombinedMonOp<Individuo> mutation(mutationA,PmutationA);
              mutation.add(mutationB, PmutationB);
              //eoPropCombinedMonOp<Individuo> mutation(mutationB,PmutationB);
              mutation.add(mutationC, PmutationC);

          //Define un objeto de encapsulacion (it contains, the crossover, the crossover rate, the mutation and the mutation rate) -> 1 line
              eoSGATransform<Individuo> encapsulacion(crossover, PcruzaT, mutation, PmutationT); //0.87

          //Define el metodo de seleccion, selecciona un individuo por cada torneo (en el parentesis se define el tama�o del torneo)
              eoDetTournamentSelect<Individuo> torneo(sizeTorneo);

          //Define un "eoSelectPerc" con el torneo como parametro por defecto (permite seleccionar el mejor individuo)
              eoSelectPerc<Individuo> seleccion(torneo);

          ////Otra estrategia de reemplazo con elitismo
              eoElitism<Individuo> reemplazo(sizeElist,false); //antes 0.6

             //Para utilizar eoElitism se define un eoDetTournamentTruncate para seleccionar los individuos para el elitismo
                  eoDetTournamentTruncate<Individuo> Trunca(sizeTorneo1);// antes 2

          //Define una poblacion de Individuos
              eoPop<Individuo> poblacion;

          //Cargar la matriz de distancias, cantidad nodos anclas y total de nodos
              Fitness.guardarDisReal(DisReal, NoAnclas, nodos, radio);

          //Cargar posiciones nodos anclas
              Fitness.guardarAnclas(vecAnclas);

          //Imprime la poblaci�n
              //poblacion.printOn(std::cout);

          //Imprime un salto de linea
              std::cout<<"Imprime un salto de linea"<< std::endl;

          //Contenedor de clases
              eoCheckPoint<Individuo> PuntoChequeo(parada);

          //Se carga el contador de generaciones al objeto eoCheckpoint para contar el numero de generaciones
              PuntoChequeo.add(*generationCounter);




          //Guardar algunas estadisticas de la poblacion
              //Muestra el mejor fitness de cada generaci�n
              //eoBestFitnessStat<Individuo> Elmejor("Mejor Fitness");
              //La media y stdev
              //eoSecondMomentStats<Individuo> SegundoStat;
              //Se agrega al eoCheckPoint
              //PuntoChequeo.add(Elmejor);
              //PuntoChequeo.add(SegundoStat);
std::cout<<"Se genera la poblacion"<< std::endl;
          /**Otra forma de cargar la poblacion**/
                      //Para la inicializaci�n del cromosoma, primero se debe definir como se generaran los genes y la semilla
                      rng.reseed(seed);
                      //Se utilizara un generador uniforme, (valor min, valor max)
                      eoUniformGenerator<double> uGen(_min, _max);
                      //Crear el inicializador para los cromosomas, llamado random
                      IndiInit random(nodos*2,uGen);

                      //Llena la poblaci�n y evalua cada cromosoma
                      for(int i=0 ; i<POP_SIZE ; i++){
                          random(cromosoma);
                          Fitness(cromosoma);
                          poblacion.push_back(cromosoma);
                      }
                      //Guarda la poblacion inicial a un archivo, para usarlo como semilla se debe agregar al inicio \section{eoPop}
                      //std::string pobla = "PopInicial.txt";
                      //std::ofstream poblacion1(pobla.c_str());
                      //poblacion.printOn(poblacion1);

          // Incializa el algoritmo genetico secuencial
              eoEasyEA<Individuo> algoritmo(PuntoChequeo, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

          //Tiempo inicial
          std::cout<<"Se inicializa el tiempo ti"<< std::endl;
              gettimeofday(&ti, NULL);

          //Corre el algoritmo en la poblacion inicializada
          //std::cout<<"Comienza la ejecucion"<< std::endl;
          std::cout<<"Comienza la ejecución del algoritmo"<< std::endl;
              algoritmo(poblacion);
          //std::cout<<"Termina la ejecucion"<< std::endl;
          //Tiempo Final
          std::cout<<"Se inicializa el tiempo tf"<< std::endl;
              gettimeofday(&tf, NULL);

              std::cout << std::endl;

          //Imprime el mejor cromosoma
              //poblacion.best_element().printOn(std::cout);
              //double valor;
              //valor = poblacion.nth_element_fitness(0);

              std::cout << std::endl;
              //std::cout <<valor<< std::endl;

          //Imprime el tiempo de ejecuci�n del algoritmo
              tiempo = (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;

              //std::cout <<"Tiempo de ejecucion en milisegundos: " << tiempo << std::endl;
              std::cout <<"Tiempo de ejecucion en segundos: " << tiempo/1000.0 << std::endl;
              //std::cout <<"Tiempo de ejecucion en minutos: " << (tiempo/1000.0)/60 << std::endl;

              //std::cout << std::endl;
          //Se grafica el error y todos los nodos
              //std::string filename="generacion";
              //graphError error(filename, setGeneracion, numberGeneration, nodos, NoAnclas, _max);

            //std::cout << std::endl;
            delete generationCounter;
            _eo.fitness(poblacion.nth_element_fitness(0));
    }
};
#endif
