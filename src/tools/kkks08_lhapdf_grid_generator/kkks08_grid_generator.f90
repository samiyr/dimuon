program main
    implicit none
    external kkks08_D

    integer :: argc, NO, NH
    character(2) :: NO_string, NH_string

    double precision :: x_list(131), Q2_list(29)
    double precision :: DH(0:10)

    integer :: i_x, i_Q2
    double precision :: x, Q2, Q

    integer :: unit
    character(len = :), allocatable :: filename

    double precision :: xg, xu, xubar, xd, xdbar, xs, xsbar, xc, xcbar, xb, xbbar

    argc = command_argument_count()
    if (argc < 2) then
        print *, "Need to specify NO and NH"
        stop
    end if

    call get_command_argument(1, NO_string)
    read(NO_string, *) NO

    call get_command_argument(2, NH_string)
    read(NH_string, *) NH
    
    x_list = (/ 1.0D-4,1.4D-4,2.0D-4,3.0D-4,4.5D-4,6.7D-4,8.0d-4,&
        1.0D-3,1.4D-3,2.0D-3,3.0D-3,4.5D-3,6.7D-3,8.0d-3,&
        1.0D-2, 1.4D-2, 2.0D-2, 3.0D-2, 4.0d-2,&
        0.046d0,0.050d0,0.054d0,0.059d0,0.064d0,0.068d0,&
        0.072d0,0.076d0,0.080d0,0.085d0,0.09d0,0.095d0,0.1d0,0.105d0,&
        0.110d0,0.117d0,0.125d0,0.130d0,0.135d0,0.140d0,0.15d0,0.16d0,&
        0.17d0,0.18d0,0.19d0,0.2d0,0.21d0,0.22d0,0.23d0,0.24d0,&
        0.246d0,0.252d0,0.260d0,0.2660d0,0.272d0,0.280d0,0.29d0,0.3d0,&
        0.31d0,0.32d0,0.33d0,0.34d0,0.35d0,0.36d0,0.37d0,&
        0.38d0,0.39d0,0.40d0,0.41d0,0.42d0,0.43d0,&
        0.44d0,0.45d0,0.46d0,0.47d0,0.48d0,0.49d0,0.5d0,&
        0.51d0,0.52d0,0.53d0,0.54d0,0.55d0,0.56d0,&
        0.57d0,0.58d0,0.59d0,0.6d0,0.61d0,0.62d0,0.63d0,&
        0.64d0,0.65d0,0.66d0,0.67d0,0.68d0,0.69d0,&
        0.70d0,0.71d0,0.72d0,0.73d0,0.74d0,0.75d0,0.76d0,&
        0.77d0,0.78d0,0.79d0,0.80d0,0.81d0,0.82d0,&
        0.83d0,0.84d0,0.85d0,0.86d0,0.87d0,0.888d0,0.89d0,&
        0.90d0,0.91d0,0.92d0,0.93d0,0.94d0,0.95d0,0.96d0,&
        0.965d0,0.97d0,0.975d0,0.98d0,0.985d0,0.99d0,0.995d0,0.9999d0&
    /)
    
    Q2_list = (/ 2.25d0, 3d0, 5d0, 7d0, 10D0, 16D0, 25D0, 40D0,&
        64D0, 81d0,&
        1.0D2, 1.5D2, 2.0D2, 2.5D2, 5.0d2,&
        1.0D3, 1.8D3, 3.2D3, 5.7D3, 8315.18d0,&
        1.0D4, 1.8D4, 3.2D4, 5.7D4, 1.0d5,&
        1.8D5, 3.2D5, 5.7D5, 1.0d6&
    /)

    DH = 0.0
    
    call system("mkdir -p LHAPDF/kkks08_"//get_set_name(NO, NH))

    filename = "LHAPDF/kkks08_"//get_set_name(NO, NH)//"/kkks08_"//create_filename(NO, NH, "_0000.dat")
    unit = 1
    open(unit, file = filename, status = 'replace')

    print *, "Writing grid file to "//filename

    write(unit, "(a)") "PdfType: central"
    write(unit, "(a)") "Format: lhagrid1"
    write(unit, "(a)") "---"

    do i_x = 1, size(x_list)
        x = x_list(i_x)
        write(unit, "(E23.16E3, ' ')", advance="no") x
    end do
    write(unit, *)
    do i_Q2 = 1, size(Q2_list)
        Q2 = Q2_list(i_Q2)
        Q = sqrt(Q2)
        write(unit, "(E23.16E3, ' ')", advance="no") sqrt(Q2)
    end do
    write(unit, *)
    write(unit, "(a)") "-5 -4 -3 -2 -1 1 2 3 4 5 21"

    do i_x = 1, size(x_list)
        do i_Q2 = 1, size(Q2_list)
            x = x_list(i_x)
            Q2 = Q2_list(i_Q2)
            Q = sqrt(Q2)
            call kkks08_D(NO, NH, x, Q2, DH)

            xg = x * DH(0)
            xu = x * DH(1)
            xubar = x * DH(2)
            xd = x * DH(3)
            xdbar = x * DH(4)
            xs = x * DH(5)
            xsbar = x * DH(6)
            xc = x * DH(7)
            xcbar = x * DH(8)
            xb = x * DH(9)
            xbbar = x * DH(10)

            write(unit, *) xbbar, xcbar, xsbar, xubar, xdbar, xd, xu, xs, xc, xb, xg
        end do
    end do

    write(unit, "(a)") "---"
    
    close(unit)

    print *, "Grid complete"

    call create_info_file(NO, NH)
contains
    subroutine create_info_file(NO, NH)
        implicit none

        integer :: NO, NH, unit
        character(len = :), allocatable :: filename, set_name
        character(1) :: tab
        
        set_name = get_set_name(NO, NH)
        filename = "LHAPDF/kkks08_"//set_name//"/kkks08_"//create_filename(NO, NH, ".info")
        unit = 2
        open(unit, file = filename, status = 'replace')

        print *, "Writing info file to "//filename

        tab = achar(9)
        
        write(unit, *) "SetDesc:", tab, tab, tab, "Original kkks08 "//set_name//" fragmentation function converted to LHAPDF grid format. Note that D^* = D^*+ + D^*- ; D^+ + D^-"
        write(unit, *) "Authors:", tab, tab, tab, "T. Kneesch, B.A. Kniehl, G. Kramer, I. Schienbein"
        write(unit, *) "Reference:", tab, tab, tab, "10.1016/j.nuclphysb.2008.02.015"
        write(unit, *) "Format:", tab, tab, tab, "lhagrid1"
        write(unit, *) "DataVersion:", tab, tab, "1"
        write(unit, *) "SetType: fragfn"
        write(unit, *) "NumMembers:", tab, tab, "1"
        write(unit, *) "Particle:", tab, tab, tab, pid(NH)
        write(unit, *) "Flavors:", tab, tab, tab, "[-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]"
        write(unit, *) "OrderQCD:", tab, tab, tab, "1"
        write(unit, *) "FlavorScheme:", tab, tab, "variable"
        write(unit, *) "NumFlavors:", tab, tab, "5"
        write(unit, *) "XMin:", tab, tab, tab, tab, "1.00e-04"
        write(unit, *) "XMax:", tab, tab, tab, tab, "9.999e-01"
        write(unit, *) "QMin:", tab, tab, tab, tab, "1.50e+00"
        write(unit, *) "QMax:", tab, tab, tab, tab, "1.00e+03 # the original interpolator issues a warning for Q > 1.00e+02"
        write(unit, *) "MUp:", tab, tab, tab, tab, "0.0"
        write(unit, *) "MDown:", tab, tab, tab, tab, "0.0"
        write(unit, *) "MStrange:", tab, tab, tab, "0.0"
        write(unit, *) "MCharm:", tab, tab, tab, "1.5"
        write(unit, *) "MBottom:", tab, tab, tab, "5.0"
        write(unit, *) "AlphaS_OrderQCD:", tab, "2"
        write(unit, *) "AlphaS_Type:", tab, tab, "analytic"
        write(unit, *) "AlphaS_Lambda4:", tab, "0.321"
        write(unit, *) "AlphaS_Lambda5:", tab, "0.221"

        close(unit)

        print *, "Info complete"
    end subroutine create_info_file
    function get_set_name(NO, NH) result(name)
        implicit none
        
        integer :: NO, NH
        
        character(len = :), allocatable :: chno, chnh, name

        if ((NO.EQ.1) .and. (NH.EQ.1)) then ! Global with mass D^*
            chno='mas'
            chnh='glob_d+st'
        elseif ((NO.EQ.1) .and. (NH.EQ.2)) then ! Global with mass D^0         
            chno='mas'
            chnh='global_d0'
        elseif ((NO.EQ.1) .and. (NH.EQ.3)) then ! Global with mass D^+
            chno='mas'
            chnh='global_d+'
        elseif ((NO.EQ.1) .and. (NH.EQ.4)) then ! BELLE fit with mass D^*
            chno='mas'
            chnh='bell_d+st'
        elseif ((NO.EQ.1) .and. (NH.EQ.5)) then ! BELLE fit with mass D^0
            chno='mas'
            chnh='belle_d0_'
        elseif ((NO.EQ.1) .and. (NH.EQ.6)) then ! BELLE fit with mass D^+
            chno='mas'
            chnh='belle_d+_'
        elseif ((NO.EQ.1) .and. (NH.EQ.7)) then ! OPAL fit with mass D^*
            chno='mas'
            chnh='opal_d+st'
        elseif ((NO.EQ.1) .and. (NH.EQ.8)) then ! OPAL fit with mass D^0
            chno='mas'
            chnh='opal_d0__'
        elseif ((NO.EQ.1) .and. (NH.EQ.9)) then ! OPAL fit with mass D^+
            chno='mas'
            chnh='opal_d+__'
        elseif ((NO.EQ.1) .and. (NH.EQ.10)) then ! ALEPH fit with mass D^*
            chno='mas'
            chnh='alep_d+st'
        elseif ((NO.EQ.1) .and. (NH.EQ.11)) then ! CLEO fit with mass D^*
            chno='mas'
            chnh='cleo_d+st'
            stop
        elseif ((NO.EQ.1) .and. (NH.EQ.12)) then ! CLEO fit with mass D^0
            chno='mas'
            chnh='cleo_d0__'
        elseif ((NO.EQ.1) .and. (NH.EQ.13)) then ! CLEO fit with mass D^+
            chno='mas'
            chnh='cleo_d+__'
            stop
        elseif ((NO.EQ.2) .and. (NH.EQ.1)) then ! Global without mass D^*
            chno='m00'
            chnh='glob_d+st'
        elseif ((NO.EQ.2) .and. (NH.EQ.2)) then ! Global without mass D^0         
            chno='m00'
            chnh='global_d0'
        elseif ((NO.EQ.2) .and. (NH.EQ.3)) then ! Global without mass D^+
            chno='m00'
            chnh='global_d+'
        elseif ((NO.EQ.2) .and. (NH.EQ.4)) then ! BELLE fit without mass D^*
            chno='m00'
            chnh='bell_d+st'
        elseif ((NO.EQ.2) .and. (NH.EQ.5)) then ! BELLE fit without mass D^0
            chno='m00'
            chnh='belle_d0_'
        elseif ((NO.EQ.2) .and. (NH.EQ.6)) then ! BELLE fit without mass D^+
            chno='m00'
            chnh='belle_d+_'
        elseif ((NO.EQ.2) .and. (NH.EQ.7)) then ! OPAL fit without mass D^*
            chno='m00'
            chnh='opal_d+st'
        elseif ((NO.EQ.2) .and. (NH.EQ.8)) then ! OPAL fit without mass D^0
            chno='m00'
            chnh='opal_d0__'
        elseif ((NO.EQ.2) .and. (NH.EQ.9)) then ! OPAL fit without mass D^+
            chno='m00'
            chnh='opal_d+__'
        elseif ((NO.EQ.2) .and. (NH.EQ.10)) then ! ALEPH fit without mass D^*
            chno='m00'
            chnh='alep_d+st'
        elseif ((NO.EQ.2) .and. (NH.EQ.11)) then ! CLEO fit without mass D^*
            chno='m00'
            chnh='cleo_d+st'
            stop
        elseif ((NO.EQ.2) .and. (NH.EQ.12)) then ! CLEO fit without mass D^0
            chno='m00'
            chnh='cleo_d0__'
            stop
        elseif ((NO.EQ.2) .and. (NH.EQ.13)) then ! CLEO fit without mass D^+
            chno='m00'
            chnh='cleo_d+__'
            stop
        else
            print*,'NO =',NO,' and ','NH =',NH,' not supported'
            stop
        end if

        name = chnh//'_'//chno
    end function get_set_name
    function create_filename(NO, NH, extension) result(filename)
        implicit none
        
        integer :: NO, NH
        
        character(*) :: extension
        character(len = :), allocatable :: filename, set_name

        set_name = get_set_name(NO, NH)
        filename = set_name//extension
    end function create_filename
    function pid(NH) result(id)
        implicit none

        integer :: NH
        character(:), allocatable :: id

        if (NH == 1 .or. NH == 4 .or. NH == 7 .or. NH == 10 .or. NH == 11) then
            id = "413" ! D^*
        elseif (NH == 2 .or. NH == 5 .or. NH == 8 .or. NH == 12) then
            id = "421" ! D^0
        elseif (NH == 3 .or. NH == 6 .or. NH == 9 .or. NH == 13) then
            id = "411" ! D^+
        end if
    end function pid
end program main
