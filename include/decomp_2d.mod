	  ÎU  ¹   k820309              13.0        þ[ T                                                                                                           
       decomp_2d.f90 DECOMP_2D       %       DECOMP_2D_INIT DECOMP_2D_FINALIZE DECOMP_INFO_INIT DECOMP_INFO_FINALIZE PARTITION DECOMP_2D_ABORT GET_DECOMP_INFO MYTYPE REAL_TYPE COMPLEX_TYPE MYTYPE_BYTES NX_GLOBAL NY_GLOBAL NZ_GLOBAL NRANK NPROC DECOMP_2D_COMM_CART_X DECOMP_2D_COMM_CART_Y DECOMP_2D_COMM_CART_Z DECOMP_INFO XSTART XEND XSIZE YSTART YEND YSIZE ZSTART ZEND ZSIZE gen@TRANSPOSE_X_TO_Y gen@TRANSPOSE_Y_TO_Z gen@TRANSPOSE_Z_TO_Y gen@TRANSPOSE_Y_TO_X gen@UPDATE_HALO gen@ALLOC_X gen@ALLOC_Y gen@ALLOC_Z                                                    
                                                             u #TRANSPOSE_X_TO_Y_REAL    #TRANSPOSE_X_TO_Y_COMPLEX 	   #         @     @X                                               #TRANSPOSE_X_TO_Y_REAL%SIZE    #TRANSPOSE_X_TO_Y_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                
 4             &                   &                   &                                                     D@@                                                
 5              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO    #         @     @X                             	                   #TRANSPOSE_X_TO_Y_COMPLEX%SIZE 
   #TRANSPOSE_X_TO_Y_COMPLEX%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                            
     SIZE               @                                 PRESENT           
 @@                                                 6             &                   &                   &                                                     D@@                                                 7              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO                                                          u #TRANSPOSE_Y_TO_Z_REAL    #TRANSPOSE_Y_TO_Z_COMPLEX    #         @     @X                                               #TRANSPOSE_Y_TO_Z_REAL%SIZE    #TRANSPOSE_Y_TO_Z_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                
 D             &                   &                   &                                                     D@@                                                
 E              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO    #         @     @X                                                #TRANSPOSE_Y_TO_Z_COMPLEX%SIZE    #TRANSPOSE_Y_TO_Z_COMPLEX%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                 F             &                   &                   &                                                     D@@                                                 G              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO                                                          u #TRANSPOSE_Z_TO_Y_REAL    #TRANSPOSE_Z_TO_Y_COMPLEX !   #         @     @X                                               #TRANSPOSE_Z_TO_Y_REAL%SIZE    #TRANSPOSE_Z_TO_Y_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                   @                                 SIZE               @                                 PRESENT           
@@@                                                
 T             &                   &                   &                                                     D@@                                                
 U              &                   &                   &                                                     
 @@                                    è             #DECOMP_INFO    #         @     @X                             !                   #TRANSPOSE_Z_TO_Y_COMPLEX%SIZE "   #TRANSPOSE_Z_TO_Y_COMPLEX%PRESENT #   #SRC $   #DST %   #OPT_DECOMP &                 @                            "     SIZE               @                            #     PRESENT           
@@@                             $                    V             &                   &                   &                                                     D@@                             %                    W              &                   &                   &                                                     
 @@                              &     è             #DECOMP_INFO                                                          u #TRANSPOSE_Y_TO_X_REAL '   #TRANSPOSE_Y_TO_X_COMPLEX -   #         @     @X                            '                   #TRANSPOSE_Y_TO_X_REAL%SIZE (   #TRANSPOSE_Y_TO_X_REAL%PRESENT )   #SRC *   #DST +   #OPT_DECOMP ,                 @                            (     SIZE               @                            )     PRESENT           
 @@                             *                   
 d             &                   &                   &                                                     D@@                             +                   
 e              &                   &                   &                                                     
 @@                              ,     è             #DECOMP_INFO    #         @     @X                             -                   #TRANSPOSE_Y_TO_X_COMPLEX%SIZE .   #TRANSPOSE_Y_TO_X_COMPLEX%PRESENT /   #SRC 0   #DST 1   #OPT_DECOMP 2                 @                            .     SIZE               @                            /     PRESENT           
 @@                             0                    f             &                   &                   &                                                     D@@                             1                    g              &                   &                   &                                                     
 @@                              2     è             #DECOMP_INFO                                                           u #UPDATE_HALO_REAL 3   #UPDATE_HALO_COMPLEX ;   #         @     @X                             3                   #UPDATE_HALO_REAL%SIZE 4   #UPDATE_HALO_REAL%PRESENT 5   #IN 6   #OUT 7   #LEVEL 8   #OPT_DECOMP 9   #OPT_GLOBAL :                 @                            4     SIZE               @                            5     PRESENT           
 @@                             6                   
 z             &                   &                   &                                                   D @@                             7                   
 {              &                   &                   &                                                     
   @                              8                      @@                              9     è              #DECOMP_INFO               @@                              :            #         @     @X                             ;                   #UPDATE_HALO_COMPLEX%SIZE <   #UPDATE_HALO_COMPLEX%PRESENT =   #IN >   #OUT ?   #LEVEL @   #OPT_DECOMP A   #OPT_GLOBAL B                 @                            <     SIZE               @                            =     PRESENT           
 @@                             >                    ~             &                   &                   &                                                   D @@                             ?                                  &                   &                   &                                                     
   @                              @                      @@                              A     è              #DECOMP_INFO               @@                              B                                                                   u #ALLOC_X_REAL C   #ALLOC_X_COMPLEX H   #         @     @X                             C                   #ALLOC_X_REAL%PRESENT D   #VAR E   #OPT_DECOMP F   #OPT_GLOBAL G                 @                            D     PRESENT         D  @                             E                   
               &                   &                   &                                                     
 @@                              F     è             #DECOMP_INFO              
 @@                              G           #         @     @X                             H                   #ALLOC_X_COMPLEX%PRESENT I   #VAR J   #OPT_DECOMP K   #OPT_GLOBAL L                 @                            I     PRESENT         D  @                             J                                  &                   &                   &                                                     
 @@                              K     è             #DECOMP_INFO              
 @@                              L                                                                  u #ALLOC_Y_REAL M   #ALLOC_Y_COMPLEX R   #         @     @X                             M                   #ALLOC_Y_REAL%PRESENT N   #VAR O   #OPT_DECOMP P   #OPT_GLOBAL Q                 @                            N     PRESENT         D  @                             O                   
               &                   &                   &                                                     
 @@                              P     è             #DECOMP_INFO              
 @@                              Q           #         @     @X                             R                   #ALLOC_Y_COMPLEX%PRESENT S   #VAR T   #OPT_DECOMP U   #OPT_GLOBAL V                 @                            S     PRESENT         D  @                             T                                  &                   &                   &                                                     
 @@                              U     è             #DECOMP_INFO              
 @@                              V                                                                  u #ALLOC_Z_REAL W   #ALLOC_Z_COMPLEX \   #         @     @X                             W                   #ALLOC_Z_REAL%PRESENT X   #VAR Y   #OPT_DECOMP Z   #OPT_GLOBAL [                 @                            X     PRESENT         D  @                             Y                   
               &                   &                   &                                                     
 @@                              Z     è             #DECOMP_INFO              
 @@                              [           #         @     @X                             \                   #ALLOC_Z_COMPLEX%PRESENT ]   #VAR ^   #OPT_DECOMP _   #OPT_GLOBAL `                 @                            ]     PRESENT         D  @                             ^                                  &                   &                   &                                                     
 @@                              _     è             #DECOMP_INFO              
 @@                              `                                                        a                                                                                                      b                                 
        L            1275070495                                             c                                 
       " L            1275072546           @@                               d                       @@                               e                       @@                               f                       @@                               g                       @@                               h                       @@                               i                       @@                               j                       @@                               k                       @@                               l                              @               @                'è                   #XST m   #XEN n   #XSZ o   #YST p   #YEN q   #YSZ r   #ZST s   #ZEN t   #ZSZ u   #X1DIST v   #Y1DIST w   #Y2DIST x   #Z2DIST y   #X1CNTS z   #Y1CNTS {   #Y2CNTS |   #Z2CNTS }   #X1DISP ~   #Y1DISP    #Y2DISP    #Z2DISP    #X1COUNT    #Y1COUNT    #Y2COUNT    #Z2COUNT    #EVEN                  $                              m                                p          p            p                                        $                              n                               p          p            p                                        $                              o                               p          p            p                                        $                              p            $                   p          p            p                                        $                              q            0                   p          p            p                                        $                              r            <                   p          p            p                                        $                              s            H                   p          p            p                                        $                              t            T                   p          p            p                                        $                              u            `              	     p          p            p                                      $                              v            p              
               &                                                       $                              w            ¸                             &                                                       $                              x                                         &                                                       $                              y            H                            &                                                       $                              z                                        &                                                       $                              {            Ø                            &                                                       $                              |                                         &                                                       $                              }            h                            &                                                       $                              ~            °                            &                                                       $                                          ø                            &                                                       $                                          @                            &                                                       $                                                                      &                                                         $                                   Ð                          $                                   Ô                          $                                   Ø                          $                                   Ü                          $                                   à                       @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                          #         @                                                      #DECOMP_2D_INIT%PRESENT    #NX    #NY    #NZ    #P_ROW    #P_COL    #PERIODIC_BC                  @                                 PRESENT           
  @@                                                   
  @@                                                   
  @@                                                   
   @                                                   
   @                                                   
 @@                                                  '   p          p            p                          #         @                                                        #         @                                                     #DECOMP_INFO_INIT%ALLOCATED    #DECOMP_INFO_INIT%MAX    #DECOMP_INFO_INIT%MOD    #NX    #NY    #NZ    #DECOMP                   @                                 ALLOCATED               @                                 MAX               @                                 MOD           
  @@                                                   
  @@                                                   
  @@                                                   
D @@                                    è              #DECOMP_INFO    #         @                                  ¡                    #DECOMP ¢             
D  @                              ¢     è              #DECOMP_INFO    #         @                                  £                    #NX ¤   #NY ¥   #NZ ¦   #PDIM §   #LSTART ¨   #LEND ©   #LSIZE ª             
   @                              ¤                     
   @                              ¥                     
   @                              ¦                     
   @                              §                    (   p          p            p                                    D  @                              ¨                    )    p          p            p                                    D  @                              ©                    *    p          p            p                                    D  @                              ª                    +    p          p            p                          #         @                                  «                    #ERRORCODE ¬   #MSG ­             
@ @@                              ¬                     
   @                             ­                    1 #         @                                   ®                    #DECOMP ¯             D  @                              ¯     è              #DECOMP_INFO                  fn#fn    À   ã  b   uapp(DECOMP_2D    £  @   J  MPI %   ã  y       gen@TRANSPOSE_X_TO_Y &   \  ­      TRANSPOSE_X_TO_Y_REAL +   	  =      TRANSPOSE_X_TO_Y_REAL%SIZE .   F  @      TRANSPOSE_X_TO_Y_REAL%PRESENT *     ¼   a   TRANSPOSE_X_TO_Y_REAL%SRC *   B  ¼   a   TRANSPOSE_X_TO_Y_REAL%DST 1   þ  Y   a   TRANSPOSE_X_TO_Y_REAL%OPT_DECOMP )   W  ³      TRANSPOSE_X_TO_Y_COMPLEX .   
  =      TRANSPOSE_X_TO_Y_COMPLEX%SIZE 1   G  @      TRANSPOSE_X_TO_Y_COMPLEX%PRESENT -     ¼   a   TRANSPOSE_X_TO_Y_COMPLEX%SRC -   C  ¼   a   TRANSPOSE_X_TO_Y_COMPLEX%DST 4   ÿ  Y   a   TRANSPOSE_X_TO_Y_COMPLEX%OPT_DECOMP %   X	  y       gen@TRANSPOSE_Y_TO_Z &   Ñ	  ­      TRANSPOSE_Y_TO_Z_REAL +   ~
  =      TRANSPOSE_Y_TO_Z_REAL%SIZE .   »
  @      TRANSPOSE_Y_TO_Z_REAL%PRESENT *   û
  ¼   a   TRANSPOSE_Y_TO_Z_REAL%SRC *   ·  ¼   a   TRANSPOSE_Y_TO_Z_REAL%DST 1   s  Y   a   TRANSPOSE_Y_TO_Z_REAL%OPT_DECOMP )   Ì  ³      TRANSPOSE_Y_TO_Z_COMPLEX .     =      TRANSPOSE_Y_TO_Z_COMPLEX%SIZE 1   ¼  @      TRANSPOSE_Y_TO_Z_COMPLEX%PRESENT -   ü  ¼   a   TRANSPOSE_Y_TO_Z_COMPLEX%SRC -   ¸  ¼   a   TRANSPOSE_Y_TO_Z_COMPLEX%DST 4   t  Y   a   TRANSPOSE_Y_TO_Z_COMPLEX%OPT_DECOMP %   Í  y       gen@TRANSPOSE_Z_TO_Y &   F  ­      TRANSPOSE_Z_TO_Y_REAL +   ó  =      TRANSPOSE_Z_TO_Y_REAL%SIZE .   0  @      TRANSPOSE_Z_TO_Y_REAL%PRESENT *   p  ¼   a   TRANSPOSE_Z_TO_Y_REAL%SRC *   ,  ¼   a   TRANSPOSE_Z_TO_Y_REAL%DST 1   è  Y   a   TRANSPOSE_Z_TO_Y_REAL%OPT_DECOMP )   A  ³      TRANSPOSE_Z_TO_Y_COMPLEX .   ô  =      TRANSPOSE_Z_TO_Y_COMPLEX%SIZE 1   1  @      TRANSPOSE_Z_TO_Y_COMPLEX%PRESENT -   q  ¼   a   TRANSPOSE_Z_TO_Y_COMPLEX%SRC -   -  ¼   a   TRANSPOSE_Z_TO_Y_COMPLEX%DST 4   é  Y   a   TRANSPOSE_Z_TO_Y_COMPLEX%OPT_DECOMP %   B  y       gen@TRANSPOSE_Y_TO_X &   »  ­      TRANSPOSE_Y_TO_X_REAL +   h  =      TRANSPOSE_Y_TO_X_REAL%SIZE .   ¥  @      TRANSPOSE_Y_TO_X_REAL%PRESENT *   å  ¼   a   TRANSPOSE_Y_TO_X_REAL%SRC *   ¡  ¼   a   TRANSPOSE_Y_TO_X_REAL%DST 1   ]  Y   a   TRANSPOSE_Y_TO_X_REAL%OPT_DECOMP )   ¶  ³      TRANSPOSE_Y_TO_X_COMPLEX .   i  =      TRANSPOSE_Y_TO_X_COMPLEX%SIZE 1   ¦  @      TRANSPOSE_Y_TO_X_COMPLEX%PRESENT -   æ  ¼   a   TRANSPOSE_Y_TO_X_COMPLEX%SRC -   ¢  ¼   a   TRANSPOSE_Y_TO_X_COMPLEX%DST 4   ^  Y   a   TRANSPOSE_Y_TO_X_COMPLEX%OPT_DECOMP     ·  o       gen@UPDATE_HALO !   &  ½      UPDATE_HALO_REAL &   ã  =      UPDATE_HALO_REAL%SIZE )      @      UPDATE_HALO_REAL%PRESENT $   `  ¼   a   UPDATE_HALO_REAL%IN %     ¼   a   UPDATE_HALO_REAL%OUT '   Ø  @   a   UPDATE_HALO_REAL%LEVEL ,      Y   a   UPDATE_HALO_REAL%OPT_DECOMP ,   q   @   a   UPDATE_HALO_REAL%OPT_GLOBAL $   ±   Ã      UPDATE_HALO_COMPLEX )   t!  =      UPDATE_HALO_COMPLEX%SIZE ,   ±!  @      UPDATE_HALO_COMPLEX%PRESENT '   ñ!  ¼   a   UPDATE_HALO_COMPLEX%IN (   ­"  ¼   a   UPDATE_HALO_COMPLEX%OUT *   i#  @   a   UPDATE_HALO_COMPLEX%LEVEL /   ©#  Y   a   UPDATE_HALO_COMPLEX%OPT_DECOMP /   $  @   a   UPDATE_HALO_COMPLEX%OPT_GLOBAL    B$  g       gen@ALLOC_X    ©$        ALLOC_X_REAL %   4%  @      ALLOC_X_REAL%PRESENT !   t%  ¼   a   ALLOC_X_REAL%VAR (   0&  Y   a   ALLOC_X_REAL%OPT_DECOMP (   &  @   a   ALLOC_X_REAL%OPT_GLOBAL     É&        ALLOC_X_COMPLEX (   W'  @      ALLOC_X_COMPLEX%PRESENT $   '  ¼   a   ALLOC_X_COMPLEX%VAR +   S(  Y   a   ALLOC_X_COMPLEX%OPT_DECOMP +   ¬(  @   a   ALLOC_X_COMPLEX%OPT_GLOBAL    ì(  g       gen@ALLOC_Y    S)        ALLOC_Y_REAL %   Þ)  @      ALLOC_Y_REAL%PRESENT !   *  ¼   a   ALLOC_Y_REAL%VAR (   Ú*  Y   a   ALLOC_Y_REAL%OPT_DECOMP (   3+  @   a   ALLOC_Y_REAL%OPT_GLOBAL     s+        ALLOC_Y_COMPLEX (   ,  @      ALLOC_Y_COMPLEX%PRESENT $   A,  ¼   a   ALLOC_Y_COMPLEX%VAR +   ý,  Y   a   ALLOC_Y_COMPLEX%OPT_DECOMP +   V-  @   a   ALLOC_Y_COMPLEX%OPT_GLOBAL    -  g       gen@ALLOC_Z    ý-        ALLOC_Z_REAL %   .  @      ALLOC_Z_REAL%PRESENT !   È.  ¼   a   ALLOC_Z_REAL%VAR (   /  Y   a   ALLOC_Z_REAL%OPT_DECOMP (   Ý/  @   a   ALLOC_Z_REAL%OPT_GLOBAL     0        ALLOC_Z_COMPLEX (   «0  @      ALLOC_Z_COMPLEX%PRESENT $   ë0  ¼   a   ALLOC_Z_COMPLEX%VAR +   §1  Y   a   ALLOC_Z_COMPLEX%OPT_DECOMP +    2  @   a   ALLOC_Z_COMPLEX%OPT_GLOBAL    @2  p       MYTYPE    °2  z       REAL_TYPE    *3  z       COMPLEX_TYPE    ¤3  @       MYTYPE_BYTES    ä3  @       NX_GLOBAL    $4  @       NY_GLOBAL    d4  @       NZ_GLOBAL    ¤4  @       NRANK    ä4  @       NPROC &   $5  @       DECOMP_2D_COMM_CART_X &   d5  @       DECOMP_2D_COMM_CART_Y &   ¤5  @       DECOMP_2D_COMM_CART_Z    ä5  o      DECOMP_INFO     S7     a   DECOMP_INFO%XST     ï7     a   DECOMP_INFO%XEN     8     a   DECOMP_INFO%XSZ     '9     a   DECOMP_INFO%YST     Ã9     a   DECOMP_INFO%YEN     _:     a   DECOMP_INFO%YSZ     û:     a   DECOMP_INFO%ZST     ;     a   DECOMP_INFO%ZEN     3<     a   DECOMP_INFO%ZSZ #   Ï<     a   DECOMP_INFO%X1DIST #   c=     a   DECOMP_INFO%Y1DIST #   ÷=     a   DECOMP_INFO%Y2DIST #   >     a   DECOMP_INFO%Z2DIST #   ?     a   DECOMP_INFO%X1CNTS #   ³?     a   DECOMP_INFO%Y1CNTS #   G@     a   DECOMP_INFO%Y2CNTS #   Û@     a   DECOMP_INFO%Z2CNTS #   oA     a   DECOMP_INFO%X1DISP #   B     a   DECOMP_INFO%Y1DISP #   B     a   DECOMP_INFO%Y2DISP #   +C     a   DECOMP_INFO%Z2DISP $   ¿C  H   a   DECOMP_INFO%X1COUNT $   D  H   a   DECOMP_INFO%Y1COUNT $   OD  H   a   DECOMP_INFO%Y2COUNT $   D  H   a   DECOMP_INFO%Z2COUNT !   ßD  H   a   DECOMP_INFO%EVEN    'E         XSTART    »E         XEND    OF         XSIZE    ãF         YSTART    wG         YEND    H         YSIZE    H         ZSTART    3I         ZEND    ÇI         ZSIZE    [J  £       DECOMP_2D_INIT '   þJ  @      DECOMP_2D_INIT%PRESENT "   >K  @   a   DECOMP_2D_INIT%NX "   ~K  @   a   DECOMP_2D_INIT%NY "   ¾K  @   a   DECOMP_2D_INIT%NZ %   þK  @   a   DECOMP_2D_INIT%P_ROW %   >L  @   a   DECOMP_2D_INIT%P_COL +   ~L     a   DECOMP_2D_INIT%PERIODIC_BC #   M  H       DECOMP_2D_FINALIZE !   ZM  À       DECOMP_INFO_INIT +   N  B      DECOMP_INFO_INIT%ALLOCATED %   \N  <      DECOMP_INFO_INIT%MAX %   N  <      DECOMP_INFO_INIT%MOD $   ÔN  @   a   DECOMP_INFO_INIT%NX $   O  @   a   DECOMP_INFO_INIT%NY $   TO  @   a   DECOMP_INFO_INIT%NZ (   O  Y   a   DECOMP_INFO_INIT%DECOMP %   íO  T       DECOMP_INFO_FINALIZE ,   AP  Y   a   DECOMP_INFO_FINALIZE%DECOMP    P         PARTITION    %Q  @   a   PARTITION%NX    eQ  @   a   PARTITION%NY    ¥Q  @   a   PARTITION%NZ    åQ     a   PARTITION%PDIM !   yR     a   PARTITION%LSTART    S     a   PARTITION%LEND     ¡S     a   PARTITION%LSIZE     5T  `       DECOMP_2D_ABORT *   T  @   a   DECOMP_2D_ABORT%ERRORCODE $   ÕT  L   a   DECOMP_2D_ABORT%MSG     !U  T       GET_DECOMP_INFO '   uU  Y   a   GET_DECOMP_INFO%DECOMP 