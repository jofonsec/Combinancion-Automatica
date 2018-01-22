#include <vector>
#include <fstream>
class escenario
{
    private:
        double DisReal[5000][5000];
        double DisRSSI[5000][5000];
        double anclas[5000];
        int nodos;
        int NoAnclas;
        std::vector<double> xt, yt;
    public:
        //Genera un escenario de forma aleatoria
        escenario(double _minimo, double _maximo, double _radio, int _nodos, int _aleatorio, int _NoAnclas); //constructor de la clase
        //Utiliza un escenario ya definido
        escenario(double _radio, int _nodos, int _NoAnclas);
        //Carga un escenario desde un archivo texto
        escenario(int _nodos, int _NoAnclas);
        void guardarNodos();
        void guardarDisReal();
        void guardarDisRSSI();
        void guardarAnclas();
        double obtenerDisRSSI(int i, int j); //no modifica el atributo receptor
        double obtenerAnclas(int i);
        double obtenerPX(int i);
        double obtenerPY(int i);
};

