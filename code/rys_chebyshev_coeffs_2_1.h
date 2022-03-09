
  #ifndef RYS_CHEBYSHEV_COEFFS_2_1_H
  #define RYS_CHEBYSHEV_COEFFS_2_1_H
  const struct {
    const int rys_order = 2;
    const int chebyshev_order = 30;
    const double x_min = 10;
    const double x_max = 20;
    const double roots_coefficients[62] = {
       1.38403875717180618491806787E-1, //   2  1  0    1
      -2.37298268253147459063767378E-2, //   2  1  1    2
       2.99551301390402915660712209E-3, //   2  1  2    3
      -4.00018729926639134538161733E-4, //   2  1  3    4
       4.97290779584712499432584850E-5, //   2  1  4    5
      -4.72600894518726495093957804E-6, //   2  1  5    6
       9.06307953787590534024302158E-8, //   2  1  6    7
       8.89525156667957597137401351E-8, //   2  1  7    8
      -2.23487470887056500629494567E-8, //   2  1  8    9
       2.32743176110094420946665192E-9, //   2  1  9   10
      2.48605343815870957412269952E-10, //   2  1 10   11
     -1.61997320872366034781354007E-10, //   2  1 11   12
      3.61230554653628205018188745E-11, //   2  1 12   13
     -4.00951783277157716606838899E-12, //   2  1 13   14
     -2.25968254311538752312681628E-13, //   2  1 14   15
      2.09768041137880274090856201E-13, //   2  1 15   16
     -4.95543278436256395132783218E-14, //   2  1 16   17
      6.00113433345148531474761769E-15, //   2  1 17   18
      1.49226610074970839753208723E-16, //   2  1 18   19
     -2.58616369222479579012637897E-16, //   2  1 19   20
      6.56769249968378187699952467E-17, //   2  1 20   21
     -8.68398823150376996558950022E-18, //   2  1 21   22
      7.14421424654740407727842304E-21, //   2  1 22   23
      3.13161335975622624518626619E-19, //   2  1 23   24
     -8.65261322203379015060983414E-20, //   2  1 24   25
      1.24768951851507116317970536E-20, //   2  1 25   26
     -2.98560911577198138547188609E-22, //   2  1 26   27
     -3.70510266314314836733015055E-22, //   2  1 27   28
      1.13001781643032522647193521E-22, //   2  1 28   29
     -1.77458816849893193313209280E-23, //   2  1 29   30
      9.34734173004786457127074200E-25, //   2  1 30   31
       4.35398702067097779315071124E-1, //   2  2  0   32
      -7.45609187269703321325709614E-2, //   2  2  1   33
       9.35619348334957949919879966E-3, //   2  2  2   34
      -1.22078914517831664094927416E-3, //   2  2  3   35
       1.39404926333065230483985721E-4, //   2  2  4   36
      -8.45186415277369389421961391E-6, //   2  2  5   37
      -1.75105395235487537543968139E-6, //   2  2  6   38
       8.26878783141651296554962383E-7, //   2  2  7   39
      -1.93728196338997374825777296E-7, //   2  2  8   40
       2.99468718213586377410766573E-8, //   2  2  9   41
      -2.23622379150088579485928815E-9, //   2  2 10   42
     -3.65192946662979048326412083E-10, //   2  2 11   43
      1.81044177482279227920186196E-10, //   2  2 12   44
     -3.97671344653359427400302984E-11, //   2  2 13   45
      5.36152257769316774224631396E-12, //   2  2 14   46
     -2.20857998436312299470912671E-13, //   2  2 15   47
     -1.11472194553408941896146775E-13, //   2  2 16   48
      3.85437598395765755549953549E-14, //   2  2 17   49
     -7.26705652062260253171023196E-15, //   2  2 18   50
      8.07852559537007435967555608E-16, //   2  2 19   51
      5.77668662630787684383602210E-18, //   2  2 20   52
     -2.75410930894592787385944949E-17, //   2  2 21   53
      7.56626613614588781709516808E-18, //   2  2 22   54
     -1.24508701715635679115335887E-18, //   2  2 23   55
      1.07254174435077081742770341E-19, //   2  2 24   56
      9.84210025176193595339402167E-21, //   2  2 25   57
     -6.20838844452999382047970611E-21, //   2  2 26   58
      1.42511019914473116876733993E-21, //   2  2 27   59
     -1.99795328054419858834471611E-22, //   2  2 28   60
      1.01590373397876899715612826E-23, //   2  2 29   61
      3.32754006096666315627359319E-24, //   2  2 30   62
  };
    const double weights_coefficients[62] = {
       2.12391195114796420538889889E-1, //   2  1  0    1
      -3.64962264839996340329427534E-2, //   2  1  1    2
       4.65466212565236379346914050E-3, //   2  1  2    3
      -6.43815361824275820635564284E-4, //   2  1  3    4
       8.85714888314313678962868673E-5, //   2  1  4    5
      -1.13158628010311998519633509E-5, //   2  1  5    6
       1.23977602072065044577781676E-6, //   2  1  6    7
      -1.05545667293513506098268202E-7, //   2  1  7    8
       7.47899427372816113281854673E-9, //   2  1  8    9
      -1.47831293325788897062214083E-9, //   2  1  9   10
      5.90561316710588142411807055E-10, //   2  1 10   11
     -1.69682349734022037076358534E-10, //   2  1 11   12
      3.26104827463396630986676193E-11, //   2  1 12   13
     -3.59572213356901350543199532E-12, //   2  1 13   14
     -1.19014991668504131169852145E-13, //   2  1 14   15
      1.55297044831824900985601867E-13, //   2  1 15   16
     -3.71591193579026793340179114E-14, //   2  1 16   17
      4.32456874183259215676533029E-15, //   2  1 17   18
      2.13341522822191624131019783E-16, //   2  1 18   19
     -2.25244726940765261397699740E-16, //   2  1 19   20
      5.54835297476203786967473166E-17, //   2  1 20   21
     -7.15594934346739159891760583E-18, //   2  1 21   22
     -6.77129918653652542791128486E-20, //   2  1 22   23
      2.87282515320708250824824602E-19, //   2  1 23   24
     -7.83524508630914163555632756E-20, //   2  1 24   25
      1.13130541410274091316930050E-20, //   2  1 25   26
     -2.54546564640137076672141070E-22, //   2  1 26   27
     -3.52020428940740482325906560E-22, //   2  1 27   28
      1.08986134271206760214063074E-22, //   2  1 28   29
     -1.76586416110966074931845268E-23, //   2  1 29   30
      1.05754660019858488493497620E-24, //   2  1 30   31
       2.15122301895915976542353740E-2, //   2  2  0   32
      -3.78452098047883917708056860E-3, //   2  2  1   33
       5.34696017176790955089732809E-4, //   2  2  2   34
      -9.83688477054156369599897261E-5, //   2  2  3   35
       2.27750259782242781645262508E-5, //   2  2  4   36
      -5.83342952483584972265020511E-6, //   2  2  5   37
       1.44017161123948919545956636E-6, //   2  2  6   38
      -3.15789073852406784864798790E-7, //   2  2  7   39
       5.86662184288146297256617715E-8, //   2  2  8   40
      -8.81715652286493324868320663E-9, //   2  2  9   41
      9.87590915776396295720202734E-10, //   2  2 10   42
     -6.71137435385321617027769516E-11, //   2  2 11   43
      1.99059834064444997649127531E-12, //   2  2 12   44
     -1.30785842639826363639226570E-12, //   2  2 13   45
      7.90926733690204956240119377E-13, //   2  2 14   46
     -2.44126252500945321353537432E-13, //   2  2 15   47
      4.84740385432497221264975210E-14, //   2  2 16   48
     -5.71220865738484516004833452E-15, //   2  2 17   49
     -4.95403322700902362794305494E-17, //   2  2 18   50
      2.06631891508923408576141122E-16, //   2  2 19   51
     -5.34468813679871379173259238E-17, //   2  2 20   52
      6.94123554327777845464705449E-18, //   2  2 21   53
      8.95368234585818498790009965E-20, //   2  2 22   54
     -2.89422701493751914446891093E-19, //   2  2 23   55
      7.85551097751815092699774379E-20, //   2  2 24   56
     -1.13315992386567540270538090E-20, //   2  2 25   57
      2.56187916518413535537492300E-22, //   2  2 26   58
      3.51879810494571628828571827E-22, //   2  2 27   59
     -1.08974463353769293649765580E-22, //   2  2 28   60
      1.76577024862233328237872332E-23, //   2  2 29   61
     -1.05747368081876170501221860E-24, //   2  2 30   62
  };
} rys_chebyshev_coeffs_2_1;
#endif /* RYS_CHEBYSHEV_COEFFS_2_1_H */

