template <class EOT>
class evalPrueba : public eoEvalFunc<EOT>{

  public :
    void operator()(EOT& _eo){
      double operadores[6];
      double Divisor = 0.0;
      double ErrorTotal = 0.0;
      std::cout << std::endl;
      std::cout << "Individuo: " << _eo << std::endl;

      for (unsigned int i= 0; i< _eo.size(); i+=7){
        operadores[i/7]=_eo[i]*64+_eo[i+1]*32+_eo[i+2]*16+_eo[i+3]*8+_eo[i+4]*4+_eo[i+5]*2+_eo[i+6]*1;
        std::cout<<operadores[i/7]<<" ";
        Divisor += operadores[i/7];
      }
      std::cout << std::endl;
      std::cout << "Probabilidad de los operadores: " << std::endl;
      for (unsigned int i=0; i< 6; i++){
        operadores[i]=operadores[i]/Divisor;
        ErrorTotal += operadores[i];
        std::cout << operadores[i] << " ";
      }

      // Finalmente se debe agregar el fitness al cromosoma
      _eo.fitness(ErrorTotal/0.8);
    }
};
