#include <eo>
#include <Individuo.h>
#include <IndiInit.h>
#include <localizacionEvalPenal.h>
#include <escenario.h>
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


template <class EOT>
class evalGAs : public eoEvalFunc<EOT>{

  public :
    void operator()(EOT& _eo){

      /**PARAMETROS DEFINIDOS**/
      //Se definen los parametros, se leen desde el parser y le asigna el valor
          unsigned seed = time(0);
          //Datos necesarios del escenario de prueba
          double _min = 0.0;
          double _max = 20.0;
          unsigned int NoAnclas = 10;
          unsigned int nodos = 100;
          double radio = 5;

          double DisReal[500][500];
          double vecAnclas[NoAnclas*2];

      //Configuracion parametros algoritmo
          unsigned int POP_SIZE = 100;
          unsigned int numberGeneration = 1000;
          unsigned int Nc = 2 ; // "Constante del operador SBX",'C',"Parametros Algoritmo").value();
          double Alpha = 0.5 ;//parser.createParam((double)(0.5), "Alpha", "Constante del operador Aritmetico",'C',"Parametros Algoritmo").value();
          float preferencia = 0.5 ;//parser.createParam((float)(0.5), "Preferencia", "Constante del operador Uniforme, define el sesgo",'C',"Parametros Algoritmo").value();
          double Pcruza = 0.87 ;//parser.createParam((double)(0.87), "Pcruza", "Probabilidad de cruzamiento SBX",'X',"Parametros Algoritmo").value();
          double epsilon = 5 ; //parser.createParam((double)(5), "Epsilon", "Rango de mutación",'F',"Parametros Algoritmo").value();
          double P_change = 0.79 ; //parser.createParam((double)(0.79), "P_cambio", "Probabilidad de que mute el gen",'E',"Parametros Algoritmo").value();
          double PmutationA = 0.33 ; //parser.createParam((double)(0.85), "Pmutacion", "Probabilidad de mutacion de la encapsulacion de SVN y Swap",'Y',"Parametros Algoritmo").value();
          double PmutationB = 0.33 ; //parser.createParam((double)(0.85), "Pmutacion1", "Probabilidad de mutacion de SVN",'Z',"Parametros Algoritmo").value();
          double PmutationC = 0.33 ; //parser.createParam((double)(0.5), "Pmutacion2", "Probabilidad de mutacion de Swap",'W',"Parametros Algoritmo").value();
          double PmutationT = 0.33 ;
          double PcruzaA = 0.33 ;
          double PcruzaB = 0.33 ;
          double PcruzaC = 0.33 ;
          double PcruzaT = 0.33 ;
          double sizeTorneo = 8 ; //parser.createParam((double)(8), "SizeTorneo", "Tamano del torneo para seleccion de individuos",'L',"Parametros Algoritmo").value();
          double sizeElist = 2 ; //parser.createParam((double)(2), "SizeElist", "Cantidad de individuos que se conservan",'B',"Parametros Algoritmo").value();
          double sizeTorneo1 = 2 ; //parser.createParam((double)(2), "SizeTorneo1", "Tamano del torneo para seleccion de individuos del elitismo",'Q',"Parametros Algoritmo").value();
      //Parametro de tiempo
          struct timeval ti, tf;
          double tiempo;

      /**CARGAR EL ESCENARIO**/
      //Escenario
          //Lee desde archivo
          escenario *pEscenario = new escenario(nodos, NoAnclas);
          //Matriz de distancia
          for (int i=0; i<nodos ; i++)
              {for (int j=0; j<nodos; j++)DisReal[i][j] = pEscenario->obtenerDisRSSI(i,j);}
          //Posicion Nodos anclas
          for (int i=0 ; i<NoAnclas*2 ; i++)vecAnclas[i] = pEscenario->obtenerAnclas(i);

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
                  eoIncrementorParam<unsigned> generationCounter("Gen.");

              //Operador especifico para el problema
              SVN mutationA(NoAnclas, numberGeneration, nodos, _min, _max, & generationCounter);

              //Mutacion incluida en EO, permite llegar mas rapido a un fitness de 600
              SWAP mutationB;

              //mutationA.setGen(& generationCounter);
              UM mutationC(bounds, epsilon, P_change);

              //Combina operadores de mutacion con su respectivo peso
              eoPropCombinedMonOp<Individuo> mutation(mutationA,PmutationA);
              mutation.add(mutationB, PmutationB);
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
              poblacion.printOn(std::cout);

          //Imprime un salto de linea
              std::cout<< std::endl;

          //Contenedor de clases
              eoCheckPoint<Individuo> PuntoChequeo(parada);

          //Se carga el contador de generaciones al objeto eoCheckpoint para contar el numero de generaciones
              PuntoChequeo.add(generationCounter);




          //Guardar algunas estadisticas de la poblacion
              //Muestra el mejor fitness de cada generaci�n
              eoBestFitnessStat<Individuo> Elmejor("Mejor Fitness");
              //La media y stdev
              eoSecondMomentStats<Individuo> SegundoStat;
              //Se agrega al eoCheckPoint
              PuntoChequeo.add(Elmejor);
              PuntoChequeo.add(SegundoStat);

          // Incializa el algoritmo genetico secuencial
              eoEasyEA<Individuo> algoritmo(PuntoChequeo, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

          //Tiempo inicial
              gettimeofday(&ti, NULL);

          //Corre el algoritmo en la poblacion inicializada
              algoritmo(poblacion);

          //Tiempo Final
              gettimeofday(&tf, NULL);

              std::cout << std::endl;

          //Imprime el mejor cromosoma
              poblacion.best_element().printOn(std::cout);

              std::cout << std::endl;
              std::cout << std::endl;

          //Imprime el tiempo de ejecuci�n del algoritmo
              tiempo = (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;

              std::cout <<"Tiempo de ejecucion en milisegundos: " << tiempo << std::endl;
              std::cout <<"Tiempo de ejecucion en segundos: " << tiempo/1000.0 << std::endl;
              std::cout <<"Tiempo de ejecucion en minutos: " << (tiempo/1000.0)/60 << std::endl;

              std::cout << std::endl;
          //Se grafica el error y todos los nodos
              //std::string filename="generacion";
              //graphError error(filename, setGeneracion, numberGeneration, nodos, NoAnclas, _max);

            std::cout << std::endl;
            //Individuo.fitness().printOn(std::cout);

    }
};
