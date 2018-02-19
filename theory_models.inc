! parameters for the control of theory models
      integer name_len, n_params
      parameter (name_len = 80)
      parameter (n_params = 10)
      type theory_param
         sequence
         character(len=name_len) model
         real, dimension(n_params) :: parameter
         logical, dimension(n_params) :: flags
      end type theory_param

      TYPE(theory_param) theory_par

      namelist /theory/theory_par
      
      common /cmn_TC/theory_par

      

      

