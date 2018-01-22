#include <fstream>
#include <string>

#include <Individuo.h>
#include <eoOp.h>

#include "gnuplot-iostream.h"

#include "eoObject.h" //Define un nombre cuando se carga o guarda un archivo state
#include "utils/eoFileMonitor.h"

//Desde
class GnuplotMonitor : public eoFileMonitor
{
    public:
        //Constructor de clases derivadas
        //El string se utiliza para definir el archivo de la clase eoFileMonitor
        GnuplotMonitor(std::string _filename,bool _activate) : eoFileMonitor(_filename, " "),
        filename(_filename), activate(_activate){}

        GnuplotMonitor(std::string _filename,bool _activate,bool _islas) : eoFileMonitor(_filename, " "),
        filename(_filename), activate(_activate), islas(_islas){}

        //Funcion virtual
        virtual ~GnuplotMonitor(){}

        virtual eoMonitor& operator()();

        //Funcion inline, inidica que cada llamado a esta función debe ser reemplazo por este codigo
        //inline void mysleep(unsigned millis) { ::usleep(millis * 1000);}

		void setGen(eoIncrementorParam<unsigned> * gen);

    private:
        eoIncrementorParam<unsigned>  * nbGenerations; /*!< numero de generaciones*/
        std::string filename;
        Gnuplot gp;
        bool activate;
        bool islas;
};
