#include <smp>
#include <eo>
#include <ga.h>
#include <evalGAs.h>
#include <evalPrueba.h>
#include <ceroOperador.h>
#include <escenario.h>
#include <sys/time.h>
typedef eoBit<eoMinimizingFitness> IndiBinario; //Es una representacion binaria donde el fitness es un double

using namespace paradiseo::smp;


int main (int argc, char* argv[]){

    //Primero se debe definir un parser que lee desde la linea de comandos o un archivo
        eoParser parser(argc, argv);
    //Se definen los parametros, se leen desde el parser y le asigna el valor
        unsigned seed = parser.createParam(unsigned(time(0)), "Semilla", "semilla de numeros aleatorios", 'S').value();
    //Configuracion parametros algoritmo
        unsigned int POP_SIZE = parser.createParam((unsigned int)(15), "PopSize", "Tamano de la poblacion",'P',"Parametros Algoritmo").value();
        unsigned int numberGeneration = parser.createParam((unsigned int)(16), "MaxGen", "Criterio de parada, Numero maximo de generaciones",'G',"Parametros Algoritmo").value();
        unsigned int pointX = parser.createParam((unsigned int)(4), "PointX", "Cantidad de puntos de cruce",'A',"Parametros Algoritmo").value();
        unsigned int countOperator = parser.createParam((unsigned int)(8), "Operadores", "Cantidad de operadores de cruce y mutacion",'C',"Parametros Algoritmo").value();
        unsigned int Alelos = parser.createParam((unsigned int)(7), "Alelos", "Cantidad de alelos por gen/operador",'D',"Parametros Algoritmo").value();
        double Pcruza = parser.createParam((double)(0.9), "Pcruza", "Probabilidad de cruzamiento SBX",'X',"Parametros Algoritmo").value();
        double Pmutation = parser.createParam((double)(0.3), "Pmutacion", "Probabilidad de mutacion de la encapsulacion de SVN y Swap",'Y',"Parametros Algoritmo").value();
        double Pmutation1 = parser.createParam((double)(0.9), "Pmutacion1", "Probabilidad de mutacion de SVN",'Z',"Parametros Algoritmo").value();
        double Pmutation2 = parser.createParam((double)(0.5), "Pmutacion2", "Probabilidad de mutacion de Swap",'W',"Parametros Algoritmo").value();
        double sizeTorneo = parser.createParam((double)(4), "SizeTorneo", "Tamano del torneo para seleccion de individuos",'L',"Parametros Algoritmo").value();
        double sizeElist = parser.createParam((double)(2), "SizeElist", "Cantidad de individuos que se conservan",'B',"Parametros Algoritmo").value();
        double sizeTorneo1 = parser.createParam((double)(4), "SizeTorneo1", "Tamano del torneo para seleccion de individuos del elitismo",'Q',"Parametros Algoritmo").value();

//Politicas de migracion
    unsigned int Intercambio1 = parser.createParam((unsigned int)(3), "Intercambio1", "Define cada cuantas generaciones migran individuos de la isla 1",'1',"Politicas de Migracion").value();
    unsigned int Intercambio2 = parser.createParam((unsigned int)(3), "Intercambio2", "Define cada cuantas generaciones migran individuos de la isla 2",'2',"Politicas de Migracion").value();
    unsigned int Intercambio3 = parser.createParam((unsigned int)(3), "Intercambio3", "Define cada cuantas generaciones migran individuos de la isla 3",'3',"Politicas de Migracion").value();
    unsigned int Intercambio4 = parser.createParam((unsigned int)(3), "Intercambio4", "Define cada cuantas generaciones migran individuos de la isla 4",'4',"Politicas de Migracion").value();
    unsigned int torneoIsla = parser.createParam((unsigned int)(4), "torneoIsla", "Tamano del torneo para seleccionar los individuos a migrar",'5',"Politicas de Migracion").value();
    unsigned int selecIsla = parser.createParam((unsigned int)(2), "selecIsla", "Cantidad de individuos seleccionados en el torneo",'6',"Politicas de Migracion").value();

//Parametros de guardado
    unsigned int setGeneracion = parser.createParam((unsigned int)(1), "setGeneracion", "Cada cuantas generaciones se guarda la poblacion",'T',"Guardar Datos").value();
    unsigned int setTime = parser.createParam((unsigned int)(0), "setTime", "Cada cuantos segundos se guarda la configuracion",'I',"Guardar Datos").value();

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

/**--------------------------------------------------------------**/

/**Partes comunes de todas las islas**/

//Generar una subclase de la clase de la funcion de evaluacion
    evalGAs<IndiBinario> Fitness(DisReal,vecAnclas);

//Es otro criterio de parada en el cual se define el minimo de generaciones y cuantas generaciones sin mejoras
    //eoSteadyFitContinue<IndiBinario> parada(10,2);

    /** CRUZAMIENTO**/
        // operador de dos puntos de cruce
          eoNPtsBitXover<IndiBinario> xover(pointX);

    /** MUTACION **/
        // Muta un bit por IndiBinario
          //eoDetBitFlip<IndiBinario> mutationOneBit;
          eoDetBitFlip<IndiBinario> mutationTwoBit(8);

        //Mutacion operador ceroOperador
          ceroOperador<IndiBinario> mutationCeroBit(countOperator);

        //Combina operadores de mutacion con su respectivo peso
          eoPropCombinedMonOp<IndiBinario> mutation(mutationTwoBit,Pmutation1);
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

                    /**Otra forma de cargar la poblacion**/
                    //Define una poblacion de Individuos
                        eoPop<IndiBinario> poblacion1;
                        eoPop<IndiBinario> poblacion2;
                        eoPop<IndiBinario> poblacion3;
                        eoPop<IndiBinario> poblacion4;
                        /** Cargar rng y poblacion desde un archivo**/
                            eoState inEstado;
                            inEstado.registerObject(rng);
                            inEstado.registerObject(poblacion1);

                        if (loadName != ""){
                            inEstado.load(loadName); //  carga la poblacion y la rng

                           if (poblacion1.size() < POP_SIZE)
                            std::cout << "WARNING, only " << poblacion1.size() << " individuals read in file " << loadName << "\nThe remaining " << POP_SIZE - poblacion1.size() << " will be randomly drawn" << std::endl;
                           if (poblacion1.size() > POP_SIZE)
                             {
                               std::cout << "WARNING, Load file contained too many individuals. Only the best will be retained" << std::endl;
                               poblacion1.resize(POP_SIZE);
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
                            poblacion1.append(POP_SIZE, random);
                            apply<IndiBinario>(Fitness, poblacion1);
                            poblacion2.append(POP_SIZE, random);
                            apply<IndiBinario>(Fitness, poblacion2);
                            poblacion3.append(POP_SIZE, random);
                            apply<IndiBinario>(Fitness, poblacion3);
                            poblacion4.append(POP_SIZE, random);
                            apply<IndiBinario>(Fitness, poblacion4);

                            //Guarda la poblacion inicial a un archivo, para usarlo como semilla se debe agregar al inicio \section{eoPop}
                            std::string pobla1 = "PopInicial1.txt";
                            std::string pobla2 = "PopInicial2.txt";
                            std::string pobla3 = "PopInicial3.txt";
                            std::string pobla4 = "PopInicial4.txt";
                            std::ofstream poblacionA(pobla1.c_str());
                            poblacion1.printOn(poblacionA);
                            std::ofstream poblacionB(pobla2.c_str());
                            poblacion2.printOn(poblacionB);
                            std::ofstream poblacionC(pobla3.c_str());
                            poblacion3.printOn(poblacionC);
                            std::ofstream poblacionD(pobla4.c_str());
                            poblacion4.printOn(poblacionD);

                        }


/**--------------------------------------------**/

/**Isla 1**/

//Define un objeto de encapsulaci�n (it contains, the crossover, the crossover rate, the mutation and the mutation rate) -> 1 line
    //eoSGATransform<IndiBinario> encapsulacion1(crossover, 0.9, mutation, 0.1); //0.87

//Criterio de parada
    eoGenContinue<IndiBinario> parada1(numberGeneration);

//Contenedor de clases
    eoCheckPoint<IndiBinario> PuntoChequeo1(parada1);

/*POLITICAS DE MIGRACION*/
//Debe enviar a los individuos cada cien generaciones
eoPeriodicContinue<IndiBinario> criterio(Intercambio1);//25
//Es un torneo entre 20 individuos
eoDetTournamentSelect<IndiBinario> selectOne1(torneoIsla); //6
//Selecciona 3 individuos usando el torneo determinista
eoSelectNumber<IndiBinario> who(selectOne1, selecIsla); //2

//La clase politicas de trabajo
MigPolicy<IndiBinario> migPolicy;
//Define las politicas especificadas
migPolicy.push_back(PolicyElement<IndiBinario>(who,criterio));

/*POLITICAS DE INTEGRACION*/
//Es una estrategia de reemplazo
eoPlusReplacement<IndiBinario> intPolicy;

//Cargar el valor de la generacion actual al operador de mutaci�n
    //Se inicializa el contador de generaciones
    eoIncrementorParam<unsigned> generationCounter1("Gen.");

/**--------------------------------------------**/


/** Guardar datos de la poblaci�n en archivos **/
    //Se carga el contador de generaciones al objeto eoCheckpoint para contar el n�mero de generaciones
    PuntoChequeo1.add(generationCounter1);
    //Genera un archivo para guardar parametros
    eoState estado1;
    //Guardar todo lo que necesites a la clase hija estado
    estado1.registerObject(poblacion1);
    //estado1.registerObject(parser);
    //Guarda el tiempo de ejecucion desde la primera generacion
    eoTimeCounter time1;
    PuntoChequeo1.add(time1);
    //Define cada cuantas generaciones se guarda la poblacion
    eoCountedStateSaver GuardarEstado1(setGeneracion,estado1,"IslaOne-Generation");
    //Siempre se debe agregar a la clase hija de eoCheckPoint para que se ejecute en cada generacion
    PuntoChequeo1.add(GuardarEstado1);

/**DATOS ESTADISTICOS**/
//Muestra el mejor fitness de cada generaci�n
    eoBestFitnessStat<IndiBinario> Elmejor1("Mejor Fitness");
    //La media (average) y stdev (standard deviation)
    eoSecondMomentStats<IndiBinario> SegundoStat1("Average & Stdev");
    //Se agrega al eoCheckPoint
    PuntoChequeo1.add(Elmejor1);
    PuntoChequeo1.add(SegundoStat1);
    // Guarda los parametros a un archivo
    eoFileMonitor fileMonitor1("Resumen1.xg", " ");
    PuntoChequeo1.add(fileMonitor1);
    fileMonitor1.add(time1);
    fileMonitor1.add(Elmejor1);           //Mejor fitness
    fileMonitor1.add(SegundoStat1);       //Media y desviacion estandar

//Inicializacion Isla 1
Island<eoEasyEA,IndiBinario> island1(poblacion1, intPolicy, migPolicy, PuntoChequeo1, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

/**--------------------------------------------**/

/**Isla 2**/

//Define un objeto de encapsulaci�n (it contains, the crossover, the crossover rate, the mutation and the mutation rate) -> 1 line
    //eoSGATransform<IndiBinario> encapsulacion2(crossover, 0.1, mutation, 0.9); //0.87

//Criterio de parada
    eoGenContinue<IndiBinario> parada2(numberGeneration);
//Contenedor de clases
    eoCheckPoint<IndiBinario> PuntoChequeo2(parada2);

/*POLITICAS DE MIGRACION*/
//Debe enviar a los individuos cada cien generaciones
eoPeriodicContinue<IndiBinario> criterio2(Intercambio2);//1000
//Es un torneo entre 20 individuos
eoDetTournamentSelect<IndiBinario> selectOne2(torneoIsla);
//Selecciona 3 individuos usando el torneo determinista
eoSelectNumber<IndiBinario> who2(selectOne2, selecIsla);


//La clase politicas de trabajo
MigPolicy<IndiBinario> migPolicy2;
//Define las politicas especificadas
migPolicy2.push_back(PolicyElement<IndiBinario>(who2,criterio2));

/*POLITICAS DE INTEGRACION*/
//Es una estrategia de reemplazo
eoPlusReplacement<IndiBinario> intPolicy2;

//Cargar el valor de la generacion actual al operador de mutaci�n
    //Se inicializa el contador de generaciones
    eoIncrementorParam<unsigned> generationCounter2("Gen.");
/**--------------------------------------------**/


/** Guardar datos de la poblaci�n en archivos **/
    //Se carga el contador de generaciones al objeto eoCheckpoint para contar el n�mero de generaciones
    PuntoChequeo2.add(generationCounter2);
    //Genera un archivo para guardar parametros
    eoState estado2;
    //Guardar todo lo que necesites a la clase hija estado
    estado2.registerObject(poblacion2);
    //Guarda el tiempo de ejecucion desde la primera generacion
    eoTimeCounter time2;
    PuntoChequeo2.add(time2);
    //Define cada cuantas generaciones se guarda la poblacion
    eoCountedStateSaver GuardarEstado2(setGeneracion,estado2,"IslaTwo-Generation");
    //Siempre se debe agregar a la clase hija de eoCheckPoint para que se ejecute en cada generacion
    PuntoChequeo2.add(GuardarEstado2);

/**DATOS ESTADISTICOS**/
//Muestra el mejor fitness de cada generaci�n
    eoBestFitnessStat<IndiBinario> Elmejor2("Mejor Fitness");
    //La media (average) y stdev (standard deviation)
    eoSecondMomentStats<IndiBinario> SegundoStat2("Average & Stdev");
    //Se agrega al eoCheckPoint
    PuntoChequeo2.add(Elmejor2);
    PuntoChequeo2.add(SegundoStat2);
    // Guarda los parametros a un archivo
    eoFileMonitor fileMonitor2("Resumen2.xg", " ");
    PuntoChequeo2.add(fileMonitor2);
    fileMonitor2.add(time2);
    fileMonitor2.add(Elmejor2);           //Mejor fitness
    fileMonitor2.add(SegundoStat2);       //Media y desviacion estandar

//Inicializacion Isla 2
Island<eoEasyEA,IndiBinario> island2(poblacion2, intPolicy2, migPolicy2, PuntoChequeo2, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

/**----------------------------------------------**/

/**Isla 3**/

//Criterio de parada
    eoGenContinue<IndiBinario> parada3(numberGeneration);
//Contenedor de clases
    eoCheckPoint<IndiBinario> PuntoChequeo3(parada3);

/*POLITICAS DE MIGRACION*/
//Debe enviar a los individuos cada X generacion
eoPeriodicContinue<IndiBinario> criterio3(Intercambio3);
//Es un torneo entre Y individuos
eoDetTournamentSelect<IndiBinario> selectOne3(torneoIsla);
//Selecciona X individuos usando el torneo determinista
eoSelectNumber<IndiBinario> who3(selectOne3,selecIsla);

//La clase politicas de trabajo
MigPolicy<IndiBinario> migPolicy3;
//Define las politicas especificadas
migPolicy3.push_back(PolicyElement<IndiBinario>(who3,criterio3));

/*POLITICAS DE INTEGRACION*/
//Es una estrategia de reemplazo
eoPlusReplacement<IndiBinario> intPolicy3;

//Cargar el valor de la generacion actual al operador de mutaci�n
    //Se inicializa el contador de generaciones
    eoIncrementorParam<unsigned> generationCounter3("Gen.");
/**--------------------------------------------**/


/** Guardar datos de la poblaci�n en archivos **/
    //Se carga el contador de generaciones al objeto eoCheckpoint para contar el n�mero de generaciones
    PuntoChequeo3.add(generationCounter3);
    //Genera un archivo para guardar parametros
    eoState estado3;
    //Guardar todo lo que necesites a la clase hija estado
    estado3.registerObject(poblacion3);
    //Guarda el tiempo de ejecucion desde la primera generacion
    eoTimeCounter time3;
    PuntoChequeo3.add(time3);
    //Define cada cuantas generaciones se guarda la poblacion
    eoCountedStateSaver GuardarEstado3(setGeneracion,estado3,"IslaThree-Generation");
    //Siempre se debe agregar a la clase hija de eoCheckPoint para que se ejecute en cada generacion
    PuntoChequeo3.add(GuardarEstado3);

/**DATOS ESTADISTICOS**/
//Muestra el mejor fitness de cada generaci�n
    eoBestFitnessStat<IndiBinario> Elmejor3("Mejor Fitness");
    //La media (average) y stdev (standard deviation)
    eoSecondMomentStats<IndiBinario> SegundoStat3("Average & Stdev");
    //Se agrega al eoCheckPoint
    PuntoChequeo3.add(Elmejor3);
    PuntoChequeo3.add(SegundoStat3);
    // Guarda los parametros a un archivo
    eoFileMonitor fileMonitor3("Resumen3.xg", " ");
    PuntoChequeo3.add(fileMonitor3);
    fileMonitor3.add(time3);
    fileMonitor3.add(Elmejor3);           //Mejor fitness
    fileMonitor3.add(SegundoStat3);       //Media y desviacion estandar

//Inicializacion Isla 3
Island<eoEasyEA,IndiBinario> island3(poblacion3, intPolicy3, migPolicy3, PuntoChequeo3, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

/**----------------------------------------------**/

/**Isla 4**/

//Define un objeto de encapsulaci�n (it contains, the crossover, the crossover rate, the mutation and the mutation rate) -> 1 line
    //eoSGATransform<IndiBinario> encapsulacion3(crossover, 0.1, mutation, 0.1); //0.87


//Criterio de parada
    eoGenContinue<IndiBinario> parada4(numberGeneration);
//Contenedor de clases
    eoCheckPoint<IndiBinario> PuntoChequeo4(parada4);

/*POLITICAS DE MIGRACION*/
//Debe enviar a los individuos cada cien generaciones
eoPeriodicContinue<IndiBinario> criterio4(Intercambio4);
//Es un torneo entre 20 individuos
eoDetTournamentSelect<IndiBinario> selectOne4(torneoIsla);
//Selecciona 3 individuos usando el torneo determinista
eoSelectNumber<IndiBinario> who4(selectOne4,selecIsla);

//La clase politicas de trabajo
MigPolicy<IndiBinario> migPolicy4;
//Define las politicas especificadas
migPolicy4.push_back(PolicyElement<IndiBinario>(who4,criterio4));

/*POLITICAS DE INTEGRACION*/
//Es una estrategia de reemplazo
eoPlusReplacement<IndiBinario> intPolicy4;

//Cargar el valor de la generacion actual al operador de mutaci�n
    //Se inicializa el contador de generaciones
    eoIncrementorParam<unsigned> generationCounter4("Gen.");
/**--------------------------------------------**/


/** Guardar datos de la poblaci�n en archivos **/
    //Se carga el contador de generaciones al objeto eoCheckpoint para contar el n�mero de generaciones
    PuntoChequeo4.add(generationCounter4);
    //Genera un archivo para guardar parametros
    eoState estado4;
    //Guardar todo lo que necesites a la clase hija estado
    estado4.registerObject(poblacion4);
    //Guarda el tiempo de ejecucion desde la primera generacion
    eoTimeCounter time4;
    PuntoChequeo4.add(time4);
    //Define cada cuantas generaciones se guarda la poblacion
    eoCountedStateSaver GuardarEstado4(setGeneracion,estado4,"IslaFour-Generation");
    //Siempre se debe agregar a la clase hija de eoCheckPoint para que se ejecute en cada generacion
    PuntoChequeo4.add(GuardarEstado4);

/**DATOS ESTADISTICOS**/
//Muestra el mejor fitness de cada generaci�n
    eoBestFitnessStat<IndiBinario> Elmejor4("Mejor Fitness");
    //La media (average) y stdev (standard deviation)
    eoSecondMomentStats<IndiBinario> SegundoStat4("Average & Stdev");
    //Se agrega al eoCheckPoint
    PuntoChequeo4.add(Elmejor4);
    PuntoChequeo4.add(SegundoStat4);
    // Guarda los parametros a un archivo
    eoFileMonitor fileMonitor4("Resumen4.xg", " ");
    PuntoChequeo4.add(fileMonitor4);
    fileMonitor4.add(time4);
    fileMonitor4.add(Elmejor4);           //Mejor fitness
    fileMonitor4.add(SegundoStat4);       //Media y desviacion estandar

//Inicializacion Isla 4
Island<eoEasyEA,IndiBinario> island4(poblacion4, intPolicy4, migPolicy4, PuntoChequeo4, Fitness, seleccion, encapsulacion, reemplazo, Trunca);

/**----------------------------------------------**/

///** Grafica **/
//    GnuplotMonitor grafica(InPut,graficaGnuplot,grafIsla);       //Grafica el fitness y la media
//    grafica.setGen(& generationCounter1); //Carga la generacion
//    PuntoChequeo1.add(grafica);
///**------------------------------------------**/


/**Topologia**/
//Define la topologia
Topology<Ring> topo;
//Inicia el algoritmo con la topologia definida
IslandModel<IndiBinario> model(topo);
/**--------------------------------------------**/

//Carga las islas
    model.add(island1);
    model.add(island2);
    model.add(island3);
    model.add(island4);

//Tiempo inicial
    gettimeofday(&ti, NULL);

//Corre el algoritmo de islas homogeneas
    model();

//Tiempo Final
    gettimeofday(&tf, NULL);

//Imprime el mejor cromosoma de cada poblacion
    std::cout << std::endl;
    poblacion1.best_element().printOn(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;
    poblacion2.best_element().printOn(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;
    poblacion3.best_element().printOn(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;
    poblacion4.best_element().printOn(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;

//Imprime el tiempo de ejecuci�n del algoritmo
    tiempo = (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;
    std::cout <<"Tiempo de ejecucion en segundos: " << tiempo/1000.0 << std::endl;
    std::cout <<"Tiempo de ejecucion en minutos: " << (tiempo/1000.0)/60 << std::endl;
    std::cout <<"Tiempo de ejecucion en milisegundos: " << tiempo << std::endl;

  std::cout << std::endl;
  return EXIT_SUCCESS;


}
