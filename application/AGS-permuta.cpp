#include <eo>
#include <ga.h>
#include <evalGAs.h>
#include <evalPrueba.h>
#include <ceroOperador.h>

#include <sys/time.h>
#include <escenario.h>
typedef eoBit<eoMinimizingFitness> IndiBinario; //Es una representacion binaria donde el fitness es un double

int main (int argc, char* argv[]){

//Primero se debe definir un parser que lee desde la linea de comandos o un archivo
    eoParser parser(argc, argv);
//Se definen los parametros, se leen desde el parser y le asigna el valor
    unsigned seed = parser.createParam(unsigned(time(0)), "Semilla", "semilla de numeros aleatorios", 'S').value();
//Configuracion parametros algoritmo
    unsigned int POP_SIZE = parser.createParam((unsigned int)(30), "PopSize", "Tamano de la poblacion",'P',"Parametros Algoritmo").value();
    unsigned int numberGeneration = parser.createParam((unsigned int)(200), "MaxGen", "Criterio de parada, Numero maximo de generaciones",'G',"Parametros Algoritmo").value();
    unsigned int pointX = parser.createParam((unsigned int)(4), "PointX", "Cantidad de puntos de cruce",'A',"Parametros Algoritmo").value();
    unsigned int countOperator = parser.createParam((unsigned int)(8), "Operadores", "Cantidad de operadores de cruce y mutacion",'C',"Parametros Algoritmo").value();
    unsigned int Alelos = parser.createParam((unsigned int)(7), "Alelos", "Cantidad de alelos por gen/operador",'D',"Parametros Algoritmo").value();
    double Pcruza = parser.createParam((double)(0.9), "Pcruza", "Probabilidad de cruzamiento SBX",'X',"Parametros Algoritmo").value();
    double Pmutation = parser.createParam((double)(0.2), "Pmutacion", "Probabilidad de mutacion de la encapsulacion de SVN y Swap",'Y',"Parametros Algoritmo").value();
    double Pmutation1 = parser.createParam((double)(0.5), "Pmutacion1", "Probabilidad de mutacion de SVN",'Z',"Parametros Algoritmo").value();
    double Pmutation2 = parser.createParam((double)(0.5), "Pmutacion2", "Probabilidad de mutacion de Swap",'W',"Parametros Algoritmo").value();
    double sizeTorneo = parser.createParam((double)(8), "SizeTorneo", "Tamano del torneo para seleccion de individuos",'L',"Parametros Algoritmo").value();
    double sizeElist = parser.createParam((double)(2), "SizeElist", "Cantidad de individuos que se conservan",'B',"Parametros Algoritmo").value();
    double sizeTorneo1 = parser.createParam((double)(2), "SizeTorneo1", "Tamano del torneo para seleccion de individuos del elitismo",'Q',"Parametros Algoritmo").value();
//Parametros de guardado
    unsigned int setGeneracion = parser.createParam((unsigned int)(20), "setGeneracion", "Cada cuantas generaciones se guarda la poblacion",'T',"Guardar Datos").value();

// El nombre del archivo status donde todos los parametros seran guardados
    std::string str_status = parser.ProgramName() + ".status"; // default value
    std::string statusName = parser.createParam(str_status, "status","Status file",'S', "Persistence" ).value();
    std::string loadName = parser.createParam(std::string(""), "Carga","Se restaura desde un archivo guardado",'L', "Persistence" ).value();

//Termina la ejecucion al consultar la ayuda
    if (parser.userNeedsHelp()){
        parser.printHelp(std::cout);
        exit(1);
          }

    if (statusName != ""){
        std::ofstream os(statusName.c_str());
          os << parser;	// and you can use that file as parameter file
          }

//Parametro de tiempo
    struct timeval ti, tf;
    double tiempo;

/**CARGAR EL ESCENARIO**/
    //Escenario
    unsigned int NoAnclas = 20;
    unsigned int nodos = 120;
    double DisReal[500][500];
    double vecAnclas[NoAnclas*2];
        //Lee desde archivo
        escenario *pEscenario = new escenario(nodos, NoAnclas);
        //Matriz de distancia
        for (int i=0; i<nodos ; i++)
            {for (int j=0; j<nodos; j++)DisReal[i][j] = pEscenario->obtenerDisRSSI(i,j);}
        //Posicion Nodos anclas
        for (int i=0 ; i<NoAnclas*2 ; i++)vecAnclas[i] = pEscenario->obtenerAnclas(i);


//Define la representacion (IndiBinario)
    //IndiBinario Individuo;

//Generar una subclase de la clase de la funcion de evaluacion
    evalGAs<IndiBinario> Fitness(DisReal,vecAnclas);

//Criterio de parada
    eoGenContinue<IndiBinario> parada(numberGeneration);

//Es otro criterio de parada en el cual se define el minimo de generaciones y cuantas generaciones sin mejoras
    //eoSteadyFitContinue<IndiBinario> parada(10,2);

/** CRUZAMIENTO**/
    // operador de dos puntos de cruce
      eoNPtsBitXover<IndiBinario> xover(pointX);

/** MUTACION **/
    // Muta un bit por IndiBinario
      eoDetBitFlip<IndiBinario> mutationOneBit;

    //Mutacion operador ceroOperador
      ceroOperador<IndiBinario> mutationCeroBit(countOperator);

    //Combina operadores de mutacion con su respectivo peso
      eoPropCombinedMonOp<IndiBinario> mutation(mutationOneBit,Pmutation1);
      mutation.add(mutationCeroBit, Pmutation2);

//Define un objeto de encapsulacion (it contains, the crossover, the crossover rate, the mutation and the mutation rate) -> 1 line
    eoSGATransform<IndiBinario> encapsulacion(xover, Pcruza, mutation, Pmutation); //0.87

/** SELECCION **/
//Seleccion para los operadores
//Define el metodo de seleccion, selecciona un IndiBinario por cada torneo (en el parentesis se define el tamano del torneo)
    eoDetTournamentSelect<IndiBinario> torneo(sizeTorneo);

//Define un "eoSelectPerc" con el torneo como parametro por defecto (representa el porcentaje)
    eoSelectPerc<IndiBinario> seleccion(torneo);

/** REEMPLAZO **/
//Reemplazo de los individuos en cada generacion
//estrategia de reemplazo con elitismo
    eoElitism<IndiBinario> reemplazo(sizeElist,false); //antes 0.6

   //Para utilizar eoElitism se define un eoDetTournamentTruncate para seleccionar los individuos para el elitismo
        eoDetTournamentTruncate<IndiBinario> Trunca(sizeTorneo1);// antes 2

//Define una poblacion de Individuos
    eoPop<IndiBinario> poblacion;

//Imprime la poblacion
    //poblacion.printOn(std::cout);

//Imprime un salto de linea
    std::cout<< std::endl;

//Contenedor de clases
    eoCheckPoint<IndiBinario> PuntoChequeo(parada);

//Se carga el contador de generaciones al objeto eoCheckpoint para contar el numero de generaciones
    //Se inicializa el contador de generaciones
      eoIncrementorParam<unsigned> generationCounter("Gen.");

    PuntoChequeo.add(generationCounter);

/** Cargar rng y poblacion desde un archivo**/
    eoState inEstado;
    inEstado.registerObject(rng);
    inEstado.registerObject(poblacion);


/** Guardar datos de la poblacion en archivos **/
    //Genera un archivo para guardar parametros

    //Guardar todo lo que necesites a la clase hija estado
    eoState outEstado;
    outEstado.registerObject(rng);
    outEstado.registerObject(poblacion);
    //outEstado.registerObject(parser); //Guarda la configuracion
    //Guarda el tiempo de ejecucion desde la primera generacion
    eoTimeCounter time;
    PuntoChequeo.add(time);
    //Define cada cuantas generaciones se guarda la poblacion
    eoCountedStateSaver GuardarEstado(setGeneracion,outEstado,"generacion");
    //Siempre se debe agregar a la clase hija de eoCheckPoint para que se ejecute en cada generacion
    PuntoChequeo.add(GuardarEstado);

//Guardar algunas estadisticas de la poblacion
    //Muestra el mejor fitness de cada generaci�n
    eoBestFitnessStat<IndiBinario> Elmejor("Mejor Fitness");
    //La media y stdev
    eoSecondMomentStats<IndiBinario> SegundoStat;
    //Se agrega al eoCheckPoint
    PuntoChequeo.add(Elmejor);
    PuntoChequeo.add(SegundoStat);
    // Guarda los parametros a un archivo
    eoFileMonitor fileMonitor("stats.xg", " ");
    PuntoChequeo.add(fileMonitor);
    fileMonitor.add(generationCounter); //Numero de generaciones
    fileMonitor.add(time);              //Tiempo total de ejecucion desde la primera generacion
    fileMonitor.add(Elmejor);           //Mejor fitness
    fileMonitor.add(SegundoStat);       //Media y desviacion estandar

/**Otra forma de cargar la poblacion**/

    if (loadName != ""){
        inEstado.load(loadName); //  carga la poblacion y la rng

       if (poblacion.size() < POP_SIZE)
        std::cout << "WARNING, only " << poblacion.size() << " individuals read in file " << loadName << "\nThe remaining " << POP_SIZE - poblacion.size() << " will be randomly drawn" << std::endl;
       if (poblacion.size() > POP_SIZE)
         {
           std::cout << "WARNING, Load file contained too many individuals. Only the best will be retained" << std::endl;
           poblacion.resize(POP_SIZE);
         }


    }
    else{
        //Para la inicializaci�n del cromosoma, primero se debe definir como se generaran los genes y la semilla
        rng.reseed(seed);
        //Se utilizara un generador uniforme
        eoUniformGenerator<bool> uGen;
        //Crear el inicializador para los cromosomas, llamado random y define un largo del vector
        eoInitFixedLength<IndiBinario> random(countOperator*Alelos,uGen);

        // Inicializa la poblacion desde el randomizer: necesita usar la funcion append
        poblacion.append(POP_SIZE, random);
        apply<IndiBinario>(Fitness, poblacion);
        //Guarda la poblacion inicial a un archivo, para usarlo como semilla se debe agregar al inicio \section{eoPop}
        std::string pobla = "PopInicial.txt";
        std::ofstream poblacion1(pobla.c_str());
        poblacion.printOn(poblacion1);
    }

// Incializa el algoritmo genetico secuencial
    eoEasyEA<IndiBinario> algoritmo(PuntoChequeo, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

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
  return EXIT_SUCCESS;


}
