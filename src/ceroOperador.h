#ifndef ceroOperador_h
#define ceroOperador_h


template<class EOT>
class ceroOperador: public eoMonOp<EOT>
{
public:
        // Constructor principal
        ceroOperador(int _countOperator):
        countOpe(_countOperator){}

        /// Descripcion de la clase
        virtual std::string className() const { return "Selecciona un gen al azar y lo deja en cero"; }

        bool operator()(EOT& _Individuo){
                unsigned aleatorio =eo::rng.random(countOpe);
                unsigned posicion = (aleatorio-1)*7;
                for(unsigned int i = posicion; i < posicion + 7; i++){
                        // Deja en cero el operador seleccionado
                        _Individuo[i]=0;
                }
                return true;
        }

private:
        int countOpe;

};

#endif
