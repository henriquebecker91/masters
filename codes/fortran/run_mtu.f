      program main

c*********************************************************************72
c
cc MAIN is the main program for hbm_mtu1.
c
c  Discussion:
c
c    Calls MTU1 or MTU2 over the filename passed as the second command
c    line arg (the first argument is 1 or 2, to specify if MTU1 or MTU2
c    is to be called).
c
c  Licensing:
c
c    Public Domain. Unlicense.
c
c  Modified:
c
c    TODO
c
c  Author:
c
c    Henrique Becker
c
      implicit none

      integer :: n    ! Quantity of items on the read instance.
      integer :: c    ! Capacity of the read instance.

      integer, allocatable :: p(:) ! Item profits of the read instance.
      integer, allocatable :: w(:) ! Item weights of the read instance.
      integer, allocatable :: x(:) ! Quantity of the items on solution.

      integer :: z    ! Solution value.
      integer :: mass ! Weight of the found optimal solution.

      integer :: jdim ! Size of MTU arrays, must be at least n+1.
      integer :: jck  ! Boolean. The instance should be checked?
      integer :: jfo  ! Boolean. Return exact result?
      integer :: jub  ! Output. Upper bound.

      ! Five arrays used internally by MTU2.
      real, allocatable :: xo(:)
      real, allocatable :: rr(:)
      integer, allocatable :: po(:)
      integer, allocatable :: wo(:)
      integer, allocatable :: pp(:)

      ! Variables needed to read the time used
      integer :: t1, t2, clock_rate, clock_max, max_sec
      real    :: ttime

      ! Variables used to parse the command-line
      character(len=255) :: mtu_version
      character(len=255) :: fname
      integer :: argc

      ! Auxiliar variables loop variables used inside main
      integer :: j
      integer :: i
      integer :: jpx
      integer :: kf

      argc = command_argument_count()
      if (argc .ne. 2) then
        write (*,*) 'usage: <binary_path> <1 or 2> <instance_filename>'
        write (*,*) 'The first parameter define if MTU1 or MTU2 will be
     *used. The second is a path to an instance in the .sukp format.'
        write (*,*) 'Number of parameters received: ', argc
        write (*,*) 'Aborting.'
        call exit(1)
      end if
      call get_command_argument(1, mtu_version)
      call get_command_argument(2, fname)

      call read_sukp_instance(fname, n, c, w, p)
      jdim = n + 1

      ! Code for printing the instance. CAUTION: knapsack_reorder
      ! changes the input of the algorithm.
c      write ( *, '(a)' ) ' '
c      write ( *, '(a)' ) '  MTU2 solves the UKP.'
c      write ( *, '(a)' ) ' '
c      write ( *, '(a,i6)' ) '  Knapsack capacity is ', c

c      call knapsack_reorder ( n, p, w )

c      write ( *, '(a)' ) ' '
c      write ( *, '(a)' ) '  Object  Profit    Mass'
c      write ( *, '(a)' ) ' '
c      do i = 1, n
c        write ( *, '(2x,i6,2x,i6,2x,i6)' ) 
c     &  i, p(i), w(i)
c      end do

      jfo = 1 ! We want an exact solution.
      jck = 1 ! We don't want input checking (to not mess the times).

      call system_clock (t1, clock_rate, clock_max)

      allocate( x(jdim) )
      allocate( xo(jdim) )
      allocate( rr(jdim) )
      allocate( po(jdim) )
      allocate( wo(jdim) )
      allocate( pp(jdim) )

      if (mtu_version .EQ. '1') then
        z = 0
        do j=1,n
          rr(j) = float(p(j))/float(w(j))
        end do
        do j=1,n
          pp(j) = j
        end do
        call sortr(n,rr,pp,jdim)
        do j=1,n
          i = pp(j)
          po(j) = p(i)
          wo(j) = w(i)
        end do
        call redu(n,po,wo,jdim,jpx,x)
        kf = 0
        j = jpx
  220   kf = kf + 1
        po(kf) = po(j)
        wo(kf) = wo(j)
        pp(kf) = pp(j)
        j = x(j)
        if (j .gt. 0) go to 220
        if (kf .gt. 1) go to 230
        xo(1) = c/wo(1)
        z = po(1)*xo(1)
        jub = z
        go to 240
  230   call mtu1(kf,po,wo,c,0.,z,xo,jdim,jub,x,rr)
  240   do j=1,n
          x(j) = 0
        end do
        do j=1,kf
          i = pp(j)
          x(i) = xo(j)
        end do
      else if (mtu_version .EQ. '2') then
        call mtu2 (n,p,w,c,z,x,jdim,jfo,jck,jub,po,wo,xo,rr,pp)
      else
        write (*,*) 'The first parameter define if MTU1 or MTU2 will be
     *used, it has to be 1 or 2.'
        write (*,*) 'Value of first argument: ', mtu_version
        write (*,*) 'Aborting.'
        call exit(1)
      end if

      call system_clock (t2, clock_rate, clock_max)
      ttime = real(t2 - t1)/real(clock_rate)
      max_sec = clock_max/clock_rate

      write (*, '(a)') ' '
      write (*, '(a)') ' Optimal Solution Breakdown:'
      write (*, '(a)') ' '
      write (*, '(a)') '  Object  Quant.  Profit    Mass'
      write (*, '(a)') ' '
      mass = 0
      do i = 1, n
        if ( x(i) .ge. 1 ) then
          mass = mass + w(i)*x(i)
          write (*, '(2x,i6,2x,i6,2x,i6,2x,i6)') i, x(i), p(i), w(i)
        end if
      end do
      write (*, '(a)') ' '

      write (*,'(a,i0)') 'opt: ', z
      write (*,'(a,i0)') 'y_opt: ', mass
      write (*,'(a,f19.9)') 'Seconds: ', ttime
      write (*,'(a,i0)') 'Max number of seconds before wrap: ', max_sec
      stop
      contains

      subroutine read_sukp_instance(fname, n, c, w, p)
c*********************************************************************72
c
cc read_sukp_instance Takes filename and reads a sUKP instance from it.
c
c  Licensing:
c
c    Public Domain. Unlicense.
c
c  Modified:
c
c    TODO
c
c  Author:
c
c    Henrique Becker
c
c  Parameters:
c
c    fname Input, Name of the file with an sUKP instance.
c    n Output, Size of the instance
c    c Output, Knapsack capacity of the instance
c    w Output, item weights
c    p Output, item profits
c
      implicit none

      character(len=*) :: fname
      integer :: n
      integer :: c
      integer, allocatable :: w(:)
      integer, allocatable :: p(:)
      integer :: tmp1, tmp2, i

      open(1,file=fname)
      read(1,*) n
      read(1,*) c
      allocate( w(n) )
      allocate( p(n) )
      do i=1,n
        read (1,*) tmp1,tmp2
        w(i) = tmp1
        p(i) = tmp2
      end do
      close(1)

      return
      end
      end program

      subroutine timestamp ( )
c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end

      subroutine knapsack_reorder ( n, p, w )
c*********************************************************************72
c
cc KNAPSACK_REORDER reorders knapsack data by "profit density".
c
c  Discussion:
c
c    This routine must be called to rearrange the data before calling
c    routines that handle a knapsack problem.
c
c    The "profit density" for object I is defined as P(I)/W(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 December 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input/output, integer P(N), the "profit" or value of each object.
c
c    Input/output, integer W(N), the "weight" or cost of each object.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(n)
      integer t
      integer w(n)
c
c  Rearrange the objects in order of "profit density".
c
      do i = 1, n
        do j = i+1, n
          if ( p(i) * w(j) < p(j) * w(i) ) then
            t    = p(i)
            p(i) = p(j)
            p(j) = t
            t    = w(i)
            w(i) = w(j)
            w(j) = t
          end if
        end do
      end do

      return
      end

      SUBROUTINE MTU2(N,P,W,C,Z,X,JDIM,JFO,JCK,JUB,PO,WO,XO,RR,PP)
C
C THIS SUBROUTINE SOLVES THE UNBOUNDED SINGLE KNAPSACK PROBLEM
C
C MAXIMIZE  Z = P(1)*X(1) + ... + P(N)*X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) .LE. C ,
C               X(J) .GE. 0 AND INTEGER  FOR J=1,...,N.
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE ENUMERATIVE ALGORITHM DESCRIBED IN
C SECTION  3.6.3 .
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDIM - 1 ;
C   2) P(J), W(J), C  POSITIVE INTEGERS;
C   3) MAX (W(J)) .LE. C .
C
C MTU2   CALLS  5  PROCEDURES: CHMTU2, KSMALL, MTU1, REDU AND SORTR.
C KSMALL CALLS  8  PROCEDURES: BLD, BLDF, BLDS1, DETNS1, DETNS2,
C                              FORWD, MPSORT AND SORT7.
C
C THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS
C ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MTU2.
C NO MACHINE-DEPENDENT CONSTANT IS USED.
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MTU2 NEEDS  8  ARRAYS ( P ,  W ,  X ,  PO ,  WO ,  XO ,  RR  AND
C                        PP ) OF LENGTH AT LEAST  JDIM .
C
C MEANING OF THE INPUT PARAMETERS:
C N    = NUMBER OF ITEM TYPES;
C P(J) = PROFIT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C W(J) = WEIGHT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C C    = CAPACITY OF THE KNAPSACK;
C JDIM = DIMENSION OF THE 8 ARRAYS;
C JFO  = 1 IF OPTIMAL SOLUTION IS REQUIRED,
C      = 0 IF APPROXIMATE SOLUTION IS REQUIRED;
C JCK  = 1 IF CHECK ON THE INPUT DATA IS DESIRED,
C      = 0 OTHERWISE.
C
C MEANING OF THE OUTPUT PARAMETERS:
C Z    = VALUE OF THE SOLUTION FOUND IF  Z .GT. 0 ,
C      = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI-
C        TION  - Z  IS VIOLATED;
C X(J) = NUMBER OF ITEMS OF TYPE  J  IN THE SOLUTION FOUND;
C JUB  = UPPER BOUND ON THE OPTIMAL SOLUTION VALUE (TO EVALUATE Z
C        WHEN JFO=0).
C
C ARRAYS PO, WO, XO, RR AND PP ARE DUMMY.
C
C ALL THE PARAMETERS BUT XO AND RR ARE INTEGER. ON RETURN OF MTU2
C ALL THE INPUT PARAMETERS ARE UNCHANGED.
C
      INTEGER P(JDIM),W(JDIM),X(JDIM),PO(JDIM),WO(JDIM),PP(JDIM),C,Z
      REAL    RR(JDIM),XO(JDIM)
      Z = 0
      IF ( JCK .EQ. 1 ) CALL CHMTU2(N,P,W,C,Z,JDIM)
      IF ( Z .LT. 0 ) RETURN
C
C HEURISTIC SOLUTION THROUGH THE CORE PROBLEM.
C
      KC = N
      IF ( N .LE. 200 ) GO TO 180
      DO 10 J=1,N
        RR(J) = FLOAT(P(J))/FLOAT(W(J))
   10 CONTINUE
      KC = N/100
      IF ( KC .LT. 100 ) KC = 100
      IF ( JFO .EQ. 0 ) KC = 100
      KK = N - KC + 1
      CALL KSMALL(N,RR,KK,(N+5)/6,XO)
      RK = RR(KK)
      IF = 0
      IL = N + 1
      DO 30 J=1,N
        RR(J) = FLOAT(P(J))/FLOAT(W(J))
        IF ( RR(J) .LT. RK ) GO TO 30
        IF ( RR(J) .EQ. RK ) GO TO 20
        IF = IF + 1
        PP(IF) = J
        GO TO 30
   20   IL = IL - 1
        PP(IL) = J
   30 CONTINUE
      IF ( IF .EQ. 0 ) GO TO 50
      CALL SORTR(IF,RR,PP,JDIM)
      DO 40 J=1,IF
        I = PP(J)
        PO(J) = P(I)
        WO(J) = W(I)
        W(I) = - W(I)
   40 CONTINUE
   50 IF1 = IF + 1
      DO 60 J=IF1,KC
        I = PP(IL)
        PO(J) = P(I)
        WO(J) = W(I)
        W(I) = - W(I)
        PP(J) = I
        IL = IL + 1
   60 CONTINUE
      JPK = PO(KC)
      JWK = WO(KC)
C REDUCTION OF THE CORE PROBLEM.
      CALL REDU(KC,PO,WO,JDIM,JPX,X)
      K = 0
      J = JPX
   70   K = K + 1
        PO(K) = PO(J)
        WO(K) = WO(J)
        PP(K) = PP(J)
        J = X(J)
      IF ( J .GT. 0 ) GO TO 70
      IF ( K .GT. 1 ) GO TO 80
      XO(1) = C/WO(1)
      IXO1 = XO(1)
      Z = PO(1)*IXO1
      JUB = Z + (C - WO(1)*IXO1)*JPK/JWK
      PO(2) = JPK
      WO(2) = JWK
      GO TO 90
C SOLUTION OF THE REDUCED CORE PROBLEM.
   80 CALL MTU1(K,PO,WO,C,RK,Z,XO,JDIM,JUB,X,RR)
   90 IF ( JFO .EQ. 0 .OR. Z .EQ. JUB ) GO TO 140
      IP1 = PO(1)
      IP2 = PO(2)
      IW1 = WO(1)
      IW2 = WO(2)
      IP3 = PO(3)
      IW3 = WO(3)
      IF ( K .GT. 2 ) GO TO 100
      IP3 = IP2
      IW3 = IW2
  100 DO 130 J=1,N
        X(J) = 0
        IF ( W(J) .GT. 0 ) GO TO 110
        W(J) = - W(J)
        GO TO 130
  110   ICR = C - W(J)
        IS1 = ICR/IW1
        IB = P(J) + IS1*IP1 + (ICR - IS1*IW1)*IP2/IW2
        IF ( IB .LE. Z ) GO TO 130
        ICRR = ICR - IS1*IW1
        IS2 = ICRR/IW2
        IP = P(J) + IS1*IP1 + IS2*IP2
        ICWS = ICRR - IS2*IW2
        IB = IP + ICWS*IP3/IW3
        ITRUNC = (IW2 - ICWS + IW1 - 1)/IW1
        LIM12 = IP + (ICWS + ITRUNC*IW1)*IP2/IW2 - ITRUNC*IP1
        IF ( LIM12 .GT. IB ) IB = LIM12
        IF ( IB .LE. Z ) GO TO 130
        DO 120 JJ=J,N
          W(JJ) = IABS(W(JJ))
  120   CONTINUE
        GO TO 180
  130 CONTINUE
      GO TO 160
  140 DO 150 J=1,N
        X(J) = 0
        W(J) = IABS(W(J))
  150 CONTINUE
  160 DO 170 J=1,K
        I = PP(J)
        X(I) = XO(J)
  170 CONTINUE
      JUB = Z
      RETURN
C
C SOLUTION THROUGH COMPLETE SORTING.
C
  180 DO 190 J=1,KC
        RR(J) = FLOAT(P(J))/FLOAT(W(J))
  190 CONTINUE
      DO 200 J=1,N
        PP(J) = J
  200 CONTINUE
      CALL SORTR(N,RR,PP,JDIM)
      DO 210 J=1,N
        I = PP(J)
        PO(J) = P(I)
        WO(J) = W(I)
  210 CONTINUE
C REDUCTION OF THE PROBLEM.
      CALL REDU(N,PO,WO,JDIM,JPX,X)
      KF = 0
      J = JPX
  220   KF = KF + 1
        PO(KF) = PO(J)
        WO(KF) = WO(J)
        PP(KF) = PP(J)
        J = X(J)
      IF ( J .GT. 0 ) GO TO 220
      IF ( KF .GT. 1 ) GO TO 230
      XO(1) = C/WO(1)
      IXO1 = XO(1)
      Z = PO(1)*IXO1
      JUB = Z
      GO TO 240
  230 CALL MTU1(KF,PO,WO,C,0.,Z,XO,JDIM,JUB,X,RR)
  240 DO 250 J=1,N
        X(J) = 0
  250 CONTINUE
      DO 260 J=1,KF
        I = PP(J)
        X(I) = XO(J)
  260 CONTINUE
      RETURN
      END
      SUBROUTINE CHMTU2(N,P,W,C,Z,JDIM)
C
C CHECK THE INPUT DATA.
C
      INTEGER P(JDIM),W(JDIM),C,Z
      IF ( N .GT. 1 .AND. N .LE. JDIM - 1 ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 DO 40 J=1,N
        IF ( P(J) .LE. 0 ) GO TO 20
        IF ( W(J) .LE. 0 ) GO TO 20
        IF ( W(J) .GT. C ) GO TO 50
   40 CONTINUE
      RETURN
   50 Z = - 3
      RETURN
      END
      SUBROUTINE KSMALL(N,S,K,N6,SS)
C SUBROUTINE TO FIND THE  K-TH  SMALLEST OF  N  ELEMENTS
C IN  O(N)  TIME.
C ENTRANCE TO KSMALL IS ACHIEVED BY USING THE STATEMENT
C        CALL KSMALL(N,S,K,(N+5)/6,SS)
C
C THE VALUES OF THE FIRST THREE PARAMETERS MUST BE DEFINED
C BY THE USER PRIOR TO CALLING KSMALL. KSMALL NEEDS ONE
C ARRAY  ( S )  OF LENGTH  N  AND ONE ARRAY ( SS )  OF
C LENGTH  (N+5)/6  .THESE ARRAYS MUST BE DIMENSIONED BY THE
C USER IN THE CALLING PROGRAM.
C
C KSMALL CALLS EIGHT SUBROUTINES: BLD, BLDF, BLDS1, DETNS1,
C DETNS2, FORWD, MPSORT AND SORT7.
C THESE SUBROUTINES ARE COMPLETELY LOCAL, I.E. THE INFORMA-
C TION THEY NEED IS PASSED THROUGH THE PARAMETER LIST.
C THE WHOLE PROGRAM IS COMPLETELY SELF CONTAINED AND COMMU-
C NICATION WITH IT IS ACHIEVED SOLELY THROUGH THE PARAMETER
C LIST OF KSMALL. NO MACHINE DEPENDENT COSTANTS ARE USED.
C THE PROGRAM IS WRITTEN IN AMERICAN NATIONAL STANDARD
C FORTRAN AND IS ACCEPTED BY THE PFORT VERIFIER.
C THE PROGRAM HAS BEEN TESTED ON A  CDC CYBER 76 , ON A  CDC
C CYBER 730  AND ON A  DIGITAL VAX 11/780 .
C
C MEANING OF THE INPUT PARAMETERS:
C N = NUMBER OF ELEMENTS.
C S = ARRAY CONTAINING THE ELEMENTS.
C K = INTEGER VALUE INDICATING THAT THE  K-TH  SMALLEST
C     ELEMENT OF  S  MUST BE FOUND ( 1 .LE. K .LE. N ) .
C
C N ,  K  AND  N6  ARE INTEGER,  S  AND  SS  ARE REAL.
C ON RETURN, THE ELEMENTS OF  S  ARE REARRANGED SO THAT
C S(K)  CONTAINS THE  K-TH  SMALLEST ELEMENT OF  S , WHILE
C S(I) .LE. S(K)  IF  I .LT. K ,  S(I) .GE. S(K)  IF
C I .GT. K .
C THE CURRENT DIMENSIONS OF WORK ARRAYS  L ,  R  AND  T
C ALLOW USE OF THE CODE FOR PRACTICALLY ANY VALUE OF  N
C ( N .LT. 98*7**5 ).
C
C IN THE FOLLOWING, THE COMMENT SECTIONS REFER TO PROCEDURE
C KSMALL DESCRIBED IN  "A HYBRID ALGORITHM FOR FINDING
C THE  K-TH  SMALLEST OF  N  ELEMENTS IN  O(N)  TIME" , BY
C M. FISCHETTI AND S. MARTELLO, ANNALS OF OPERATIONAL
C RESEARCH 13, 1988.
C
      DIMENSION S(N),SS(N6)
      INTEGER L(6),R(6),T(6)
      L(1) = 1
      R(1) = N
      T(1) = K
      LEV = 1
C
C STATEMENTS  1 - 10 .
C
   10 IF ( LEV .GT. 1 )  GO TO 20
      CALL FORWD(N,S,N6,SS,LEV,L,R,T,1,V,JFLAG)
      GO TO 30
   20 CALL FORWD(N6,SS,N6,SS,LEV,L,R,T,R(LEV)+1,V,JFLAG)
   30 IF ( JFLAG .EQ. 0 ) GO TO 20
   40 LEV = LEV - 1
      IF ( LEV .EQ. 0 ) RETURN
      IL = L(LEV)
      IR = R(LEV)
      IT = T(LEV)
      NN = IR - IL + 1
      NN7 = NN/7
      II = IL + 7*(NN7 - 1)
      NXT = 1
      IF ( LEV .GT. 1 ) NXT = IR + 1
C
C STATEMENTS  11 - 13 .
C
C COMPUTE  NS1 = CARDINALITY OF SET  A1 .
      IF ( LEV .GT. 1 )  GO TO 50
      CALL DETNS1(N,S,N6,SS,IL,IR,II,NXT,V,NS1,NLR)
      GO TO 60
   50 CALL DETNS1(N6,SS,N6,SS,IL,IR,II,NXT,V,NS1,NLR)
   60 IF ( NS1 .LT. IT ) GO TO 90
C EXPLICITLY DETERMINE SET  A1 .
      IF ( LEV .GT. 1 ) GO TO 70
      CALL BLDF(N,S,N6,SS,IL,IR,NLR)
      GO TO 80
   70 CALL BLDS1(N6,SS,IL,II,NXT,NLR)
   80 R(LEV) = IL + NS1 - 1
      GO TO 10
   90 IF ( NS1 .LT. 11*NN/70 ) GO TO 110
C EXPLICITLY DETERMINE SET  A - A1 .
      T(LEV) = IT - NS1
      IF ( LEV .GT. 1 ) GO TO 100
      CALL BLDF(N,S,N6,SS,IL,IR,NLR)
      L(LEV) = IL + NS1
      GO TO 10
  100 CALL BLD(N6,SS,IL,IR,II,NXT,NLR)
      R(LEV) = IL + ( NN - NS1 ) - 1
      GO TO 10
C
C STATEMENTS  14 - 16 .
C
C COMPUTE  NS2 = CARDINALITY OF SET  A2 .
  110 IF ( LEV .GT. 1 ) GO TO 120
      CALL DETNS2(N,S,N6,SS,IL,IR,II,NXT,V,NS2,NLR)
      GO TO 130
  120 CALL DETNS2(N6,SS,N6,SS,IL,IR,II,NXT,V,NS2,NLR)
  130 NS12 = NS1 + NS2
      IF ( NS12 .LT. IT ) GO TO 160
      IF ( LEV .GT. 1 ) GO TO 40
      CALL BLDF(N,S,N6,SS,IL,IR,NLR)
      DO 140 I=1,NS12
        IF ( S(I) .EQ. V ) GO TO 150
  140 CONTINUE
  150 S(I) = S(K)
      S(K) = V
      RETURN
C EXPLICITLY DETERMINE SET  A - A1 - A2 .
  160 T(LEV) = IT - NS12
      IF ( LEV .GT. 1 ) GO TO 170
      CALL BLDF(N,S,N6,SS,IL,IR,NLR)
      L(LEV) = IL + NS12
      GO TO 10
  170 CALL BLD(N6,SS,IL,IR,II,NXT,NLR)
      R(LEV) = IL + ( NN - NS12 ) - 1
      GO TO 10
      END
      SUBROUTINE BLD(N6,SS,IL,IR,II,NXT,NLR)
C SUBROUTINE TO EXPLICITLY DETERMINE SET  A - A1  OR SET
C A - A1 - A2  WHEN LEV .GT. 1 .
      DIMENSION SS(N6)
      J = IL - 1
      NEXT = NXT
      DO 30 I=IL,II,7
        IK = SS(NEXT) + 1.
        I6 = I + 6
        IF ( IK .GT. I6 ) GO TO 20
        DO 10 IJ=IK,I6
          J = J + 1
          SS(J) = SS(IJ)
   10   CONTINUE
   20   NEXT = NEXT + 1
   30 CONTINUE
      IRR = II + 7 + NLR
      IF ( IRR .GT. IR ) RETURN
      DO 40 IJ=IRR,IR
        J = J + 1
        SS(J) = SS(IJ)
   40 CONTINUE
      RETURN
      END
      SUBROUTINE BLDF(N,S,N6,SS,IL,IR,NLR)
C SUBROUTINE TO EXPLICITLY DETERMINE SET  A1  OR SET  A - A1
C OR SET  A - A1 - A2  WHEN  LEV .EQ. 1 .
      DIMENSION S(N),SS(N6)
      NN7 = (IR - IL + 1)/7
      ICS = 1
      ICD = NN7 + 1
      IUS = IL + 6
      IPS = SS(1) + 1.
      IUD = IL + 7*NN7
      IPD = IUD + NLR - 1
   10 IF ( IPS .LE. IUS ) GO TO 20
      ICS = ICS + 1
      IF ( ICS .EQ. ICD ) RETURN
      IUS = IUS + 7
      IPS = SS(ICS) + 1.
      GO TO 10
   20 IF ( IPD .GE. IUD ) GO TO 30
      ICD = ICD - 1
      IF ( ICD .EQ. ICS ) RETURN
      IUD = IUD - 7
      IPD = SS(ICD)
      GO TO 20
   30 AP = S(IPS)
      S(IPS) = S(IPD)
      S(IPD) = AP
      IPS = IPS + 1
      IPD = IPD - 1
      GO TO 10
      END
      SUBROUTINE BLDS1(N6,SS,IL,II,NXT,NLR)
C SUBROUTINE TO EXPLICITLY DETERMINE SET  A1  WHEN
C LEV .GT. 1 .
      DIMENSION SS(N6)
      J = IL - 1
      NEXT = NXT
      DO 30 I=IL,II,7
        IK = SS(NEXT)
        IF ( IK .LT. I ) GO TO 20
        DO 10 IJ=I,IK
          J = J + 1
          SS(J) = SS(IJ)
   10   CONTINUE
   20   NEXT = NEXT + 1
   30 CONTINUE
      IK = II + 7
      IRR = IK + NLR - 1
      IF ( IK .GT. IRR ) RETURN
      DO 40 IJ=IK,IRR
        J = J + 1
        SS(J) = SS(IJ)
   40 CONTINUE
      RETURN
      END
      SUBROUTINE DETNS1(NA,A,N6,SS,IL,IR,II,NXT,V,NS1,NLR)
C SUBROUTINE TO COMPUTE THE CARDINALITY OF SET  A1 .
      DIMENSION A(NA),SS(N6)
      NS1 = 0
      NEXT = NXT
      DO 90 I=IL,II,7
        IF ( A(I+3) .LT. V ) GO TO 40
        IF ( A(I+1) .LT. V ) GO TO 20
        IF ( A(I) .LT. V ) GO TO 10
        SS(NEXT) = I - 1
        GO TO 80
   10   NS1 = NS1 + 1
        SS(NEXT) = I
        GO TO 80
   20   IF ( A(I+2) .LT. V ) GO TO 30
        NS1 = NS1 + 2
        SS(NEXT) = I + 1
        GO TO 80
   30   NS1 = NS1 + 3
        SS(NEXT) = I + 2
        GO TO 80
   40   IF ( A(I+5) .LT. V ) GO TO 60
        IF ( A(I+4) .LT. V ) GO TO 50
        NS1 = NS1 + 4
        SS(NEXT) = I + 3
        GO TO 80
   50   NS1 = NS1 + 5
        SS(NEXT) = I + 4
        GO TO 80
   60   IF ( A(I+6) .LT. V ) GO TO 70
        NS1 = NS1 + 6
        SS(NEXT) = I + 5
        GO TO 80
   70   NS1 = NS1 + 7
        SS(NEXT) = I + 6
   80   NEXT = NEXT + 1
   90 CONTINUE
      NLR = 0
      IRR = II + 7
      IF ( IRR .GT. IR ) RETURN
      DO 100 I=IRR,IR
        IF ( A(I) .GE. V ) GO TO 110
        NLR = NLR + 1
  100 CONTINUE
  110 NS1 = NS1 + NLR
      RETURN
      END
      SUBROUTINE DETNS2(NA,A,N6,SS,IL,IR,II,NXT,V,NS2,NLR)
C SUBROUTINE TO COMPUTE THE CARDINALITY OF SET  A2 .
      DIMENSION A(NA),SS(N6)
      NS2 = 0
      NEXT = NXT
      DO 50 I=IL,II,7
        IS = SS(NEXT) + 1.
        I6 = I + 6
        IF ( IS .LE. I6 ) GO TO 10
        SS(NEXT) = I6
        GO TO 40
   10   DO 20 J=IS,I6
          IF ( A(J) .GT. V ) GO TO 30
   20   CONTINUE
        NS2 = NS2 + I6 - IS + 1
        SS(NEXT) = I6
        GO TO 40
   30   NS2 = NS2 + J - IS
        SS(NEXT) = J - 1
   40   NEXT = NEXT + 1
   50 CONTINUE
      IRR = II + 7 + NLR
      IF ( IRR .GT. IR ) RETURN
      NER = 0
      DO 60 I=IRR,IR
        IF ( A(I) .GT. V ) GO TO 70
        NER = NER + 1
   60 CONTINUE
   70 NS2 = NS2 + NER
      NLR = NLR + NER
      RETURN
      END
      SUBROUTINE FORWD(NA,A,N6,SS,LEV,L,R,T,NXT,V,JFLAG)
C SUBROUTINE TO PERFORM STATEMENTS  1 - 9 .
      DIMENSION A(NA),SS(N6)
      INTEGER L(6),R(6),T(6)
      IL = L(LEV)
      IR = R(LEV)
      IT = T(LEV)
      NN = IR - IL + 1
      ITARG = NN
C STATEMENT  1 .
      IF ( NN .GT. 97 ) GO TO 20
   10 CALL MPSORT(NA,A,IL,IR,IT,V)
      JFLAG = 1
      RETURN
C STATEMENT  2 .
   20 ITARG = 59*ITARG/70
      ILC = IL
      IMED = (IL + IR)/2
      P = A(IMED)
      IF ( A(IL) .LE. P ) GO TO 30
      A(IMED) = A(IL)
      A(IL) = P
      P = A(IMED)
   30 IRC = IR
      IF ( A(IR) .GE. P ) GO TO 50
      A(IMED) = A(IR)
      A(IR) = P
      P = A(IMED)
      IF ( A(IL) .LE. P ) GO TO 50
      A(IMED) = A(IL)
      A(IL) = P
      P = A(IMED)
      GO TO 50
   40 A(IRC) = A(ILC)
      A(ILC) = AUX
   50 IRC = IRC - 1
      IF ( A(IRC) .GT. P ) GO TO 50
      AUX = A(IRC)
   60 ILC = ILC + 1
      IF ( A(ILC) .LT. P ) GO TO 60
      IF ( ILC .LE. IRC ) GO TO 40
C STATEMENT  3 .
      IF ( IT .GT. ILC - IL ) GO TO 70
      IR = ILC - 1
      GO TO 80
C STATEMENT  4 .
   70 IT = IT - ( ILC - IL )
      IL = ILC
   80 NN = IR - IL + 1
      IF ( NN .LE. 97 ) GO TO 10
C STATEMENT  5 .
      IF ( NN .LE. ITARG ) GO TO 20
C STATEMENT  6 .
      JFLAG = 0
      L(LEV) = IL
      R(LEV) = IR
      T(LEV) = IT
      LEV = LEV + 1
      NN7 = NN/7
      II = IL + 7*(NN7 - 1)
      NEXT = NXT
      L(LEV) = NXT
      DO 90 I=IL,II,7
        CALL SORT7(NA,A,I)
        SS(NEXT) = A(I+3)
        NEXT = NEXT + 1
   90 CONTINUE
      I1 = II + 7
      I2 = IR - 1
  100 IF ( I1 .GT. I2 ) GO TO 120
      DO 110 I=I1,I2
        IF ( A(I) .LE. A(I+1) ) GO TO 110
        AP = A(I)
        A(I) = A(I+1)
        A(I+1) = AP
  110 CONTINUE
      I2 = I2 - 1
      GO TO 100
  120 R(LEV) = NEXT - 1
C STATEMENT  7 .
      IT1 = IT
      N1 = (11*NN + 279)/280
      IT = IT/7
      IF ( N1 .GT. IT ) IT = N1
      N2 = NN7 - N1 + 1
      IF ( N2 .LT. IT ) IT = N2
C STATEMENT  8 .
      IF ( IT .LE. (IT1 + 3)/4 ) GO TO 130
      T(LEV) = (IT1 + 3)/4
      RETURN
C STATEMENT  9 .
  130 ITT = NN7 - (NN - IT1 + 4)/4 + 1
      IF ( IT .LT. ITT ) IT = ITT
      T(LEV) = IT
      RETURN
      END
      SUBROUTINE MPSORT(NA,A,I1,I2,IT,V)
C SUBROUTINE TO REARRANGE THE ARRAY SEGMENT  A(I1:I2) SO
C THAT  A(IT+I1-1)  CONTAINS THE  IT-TH  SMALLEST ELEMENT.
      DIMENSION A(NA)
      I = I1
      J = I2
      ITA = IT + I1 - 1
   10 IF ( I .LT. J ) GO TO 20
      V = A(ITA)
      RETURN
   20 IRR = I
      IJ = (I + J)/2
      AP = A(IJ)
      IF ( A(I) .LE. AP ) GO TO 30
      A(IJ) = A(I)
      A(I) = AP
      AP = A(IJ)
   30 ILL = J
      IF ( A(J) .GE. AP ) GO TO 50
      A(IJ) = A(J)
      A(J) = AP
      AP = A(IJ)
      IF ( A(I) .LE. AP ) GO TO 50
      A(IJ) = A(I)
      A(I) = AP
      AP = A(IJ)
      GO TO 50
   40 A(ILL) = A(IRR)
      A(IRR) = AUX
   50 ILL = ILL - 1
      IF ( A(ILL) .GT. AP ) GO TO 50
      AUX = A(ILL)
   60 IRR = IRR + 1
      IF ( A(IRR) .LT. AP ) GO TO 60
      IF ( IRR .LE. ILL ) GO TO 40
      IF ( ITA .LT. IRR ) GO TO 70
      I = IRR
      GO TO 80
   70 J = ILL
   80 IF ( J - I .GT. 10 ) GO TO 20
      IF ( I .EQ. I1 ) GO TO 10
      I = I - 1
   90 I = I + 1
      IF ( I .NE. J ) GO TO 100
      V = A(ITA)
      RETURN
  100 AP = A(I+1)
      IF ( A(I) .LE. AP ) GO TO 90
      IRR = I
  110 A(IRR+1) = A(IRR)
      IRR = IRR - 1
      IF ( AP .LT. A(IRR) ) GO TO 110
      A(IRR+1) = AP
      GO TO 90
      END
      SUBROUTINE SORT7(NA,A,I)
C SUBROUTINE TO SORT IN INCREASING ORDER THE ELEMENTS FROM
C A(I)  TO  A(I+6)  OF  A  BY PERFORMING AT MOST  13  TESTS.
      DIMENSION A(NA)
      I1 = I + 1
      I2 = I + 2
      I3 = I + 3
      I4 = I + 4
      I5 = I + 5
      I6 = I + 6
      IF ( A(I) .GT. A(I1) ) GO TO 10
      A1 = A(I)
      A2 = A(I1)
      GO TO 20
   10 A1 = A(I1)
      A2 = A(I)
   20 IF ( A(I2) .GT. A(I3) ) GO TO 30
      A3 = A(I2)
      A4 = A(I3)
      GO TO 40
   30 A3 = A(I3)
      A4 = A(I2)
   40 IF ( A1 .GT. A3 ) GO TO 50
      A5 = A2
      A2 = A3
      A3 = A4
      GO TO 60
   50 A5 = A4
      AUX = A3
      A3 = A2
      A2 = A1
      A1 = AUX
   60 A4 = A(I4)
      IF ( A4 .GE. A2 ) GO TO 80
      IF ( A4 .GE. A1 ) GO TO 70
      A4 = A3
      A3 = A2
      A2 = A1
      A1 = A(I4)
      GO TO 90
   70 A4 = A3
      A3 = A2
      A2 = A(I4)
      GO TO 90
   80 IF ( A4 .GE. A3 ) GO TO 90
      A4 = A3
      A3 = A(I4)
   90 A(I) = A1
      IF ( A5 .GT. A3 ) GO TO 110
      A(I3) = A3
      A(I4) = A4
      IF ( A5 .GT. A2 ) GO TO 100
      A(I1) = A5
      A(I2) = A2
      GO TO 130
  100 A(I1) = A2
      A(I2) = A5
      GO TO 130
  110 A(I1) = A2
      A(I2) = A3
      IF ( A5 .GT. A4 ) GO TO 120
      A(I3) = A5
      A(I4) = A4
      GO TO 130
  120 A(I3) = A4
      A(I4) = A5
  130 A5 = A(I5)
      IF ( A5 .LT. A(I2) ) GO TO 150
      IF ( A5 .LT. A(I3) ) GO TO 140
      IF ( A5 .GE. A(I4) ) GO TO 180
      A(I5) = A(I4)
      A(I4) = A5
      GO TO 180
  140 A(I5) = A(I4)
      A(I4) = A(I3)
      A(I3) = A5
      GO TO 180
  150 A(I5) = A(I4)
      A(I4) = A(I3)
      A(I3) = A(I2)
      IF ( A5 .LT. A(I1) ) GO TO 160
      A(I2) = A5
      GO TO 180
  160 A(I2) = A(I1)
      IF ( A5 .LT. A(I) ) GO TO 170
      A(I1) = A5
      GO TO 180
  170 A(I1) = A(I)
      A(I) = A5
  180 A6 = A(I6)
      IF ( A6 .LT. A(I3) ) GO TO 200
      IF ( A6 .LT. A(I4) ) GO TO 190
      IF ( A6 .GE. A(I5) ) RETURN
      A(I6) = A(I5)
      A(I5) = A6
      RETURN
  190 A(I6) = A(I5)
      A(I5) = A(I4)
      A(I4) = A6
      RETURN
  200 A(I6) = A(I5)
      A(I5) = A(I4)
      A(I4) = A(I3)
      IF ( A6 .LT. A(I1) ) GO TO 220
      IF ( A6 .LT. A(I2) ) GO TO 210
      A(I3) = A6
      RETURN
  210 A(I3) = A(I2)
      A(I2) = A6
      RETURN
  220 A(I3) = A(I2)
      A(I2) = A(I1)
      IF ( A6 .LT. A(I) ) GO TO 230
      A(I1) = A6
      RETURN
  230 A(I1) = A(I)
      A(I) = A6
      RETURN
      END
      SUBROUTINE MTU1(N,P,W,C,RN,Z,X,JDIM,JUB,XX,MIN)
C
C THIS SUBROUTINE SOLVES THE UNBOUNDED SINGLE KNAPSACK PROBLEM
C
C MAXIMIZE  Z = P(1)*X(1) + ... + P(N)*X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) .LE. C ,
C               X(J) .GE. 0 AND INTEGER  FOR J=1,...,N.
C
C THE PROGRAM IS BASED ON THE BRANCH-AND-BOUND ALGORITHM PRESENTED IN
C  S. MARTELLO, P. TOTH, "BRANCH AND BOUND ALGORITHMS FOR THE SOLUTION
C  OF THE GENERAL UNIDIMENSIONAL KNAPSACK PROBLEM", IN M. ROUBENS, ED.,
C  "ADVANCES IN OPERATIONS RESEARCH", NORTH HOLLAND, 1977.
C
      INTEGER P(JDIM),W(JDIM),C,Z
      INTEGER XX(JDIM),CWS,CWF,DIFF,R,S,S1,S2,T,PROFIT,PS
      REAL    X(JDIM),MIN(JDIM)
C
C STEP 1.
C
      CWF = C
      S1 = C/W(1)
      S2 = (C - S1*W(1))/W(2)
      IP = S1*P(1) + S2*P(2)
      CWS = C - S1*W(1) - S2*W(2)
      IF ( CWS .NE. 0 ) GO TO 20
      Z = IP
      JUB = Z
      X(1) = S1
      X(2) = S2
      IF ( N .EQ. 2 ) RETURN
      DO 10 J=3,N
        X(J) = 0
   10 CONTINUE
      RETURN
   20 MINK = C + 1
      MIN(N) = MINK
      DO 30 J=2,N
        K = N + 2 - J
        IF ( W(K) .LT. MINK ) MINK = W(K)
        MIN(K-1) = MINK
   30 CONTINUE
      W(N+1) = C + 1
      P(N+1) = 0
      LIM = IP + CWS*P(3)/W(3)
      IF ( N .EQ. 2 ) LIM = FLOAT(IP) + FLOAT(CWS)*RN
      ITRUNC = (W(2) - CWS + W(1) - 1)/W(1)
      LIM12 = IP + (CWS + ITRUNC*W(1))*P(2)/W(2) - ITRUNC*P(1)
      IF ( LIM12 .GT. LIM ) LIM = LIM12
      JUB = LIM
      Z = 0
      XX(1) = S1
      XX(2) = S2
      IF ( N .EQ. 2 ) GO TO 50
      DO 40 J=3,N
        XX(J) = 0
   40 CONTINUE
   50 PROFIT = IP
      C = CWS
      II = 2
      GO TO 110
C
C STEP 2.
C
   60 S = C/W(II)
      IF ( S .GT. 0 ) GO TO 70
      IF ( Z .GE. PROFIT + C*P(II+1)/W(II+1) ) GO TO 120
      II = II + 1
      GO TO 60
C
C STEP 3.
C
   70 PS = PROFIT + S*P(II)
      CWS = C - S*W(II)
      IF ( ( CWS .EQ. 0 ) .OR. ( II .EQ. N ) ) GO TO 80
      IF ( Z .GE. PS + CWS*P(II+1)/W(II+1) ) GO TO 150
      C = CWS
      PROFIT = PS
      XX(II) = S
      GO TO 110
   80 IF ( Z .GE. PS ) GO TO 150
      Z = PS
      II1 = II - 1
      DO 90 J=1,II1
        X(J) = XX(J)
   90 CONTINUE
      X(II) = S
      II1 = II + 1
      DO 100 J=II1,N
        X(J) = 0
  100 CONTINUE
      IF ( Z .NE. LIM ) GO TO 150
      C = CWF
      RETURN
C
C STEP 4.
C
  110 II = II + 1
      IF ( C .GE. INT(MIN(II-1)) ) GO TO 60
C
C STEP 5.
C
  120 IF ( Z .GE. PROFIT ) GO TO 140
      Z = PROFIT
      DO 130 J=1,N
        X(J) = XX(J)
  130 CONTINUE
      IF ( Z .NE. LIM ) GO TO 140
      C = CWF
      RETURN
  140 IF ( XX(N) .EQ. 0 ) GO TO 150
      C = C + XX(N)*W(N)
      PROFIT = PROFIT - XX(N)*P(N)
      XX(N) = 0
C
C STEP 6.
C
  150 IB = II - 1
      DO 160 J=1,IB
        KK = II - J
        IF ( XX(KK) .GT. 0 ) GO TO 170
  160 CONTINUE
      C = CWF
      RETURN
  170 R = C
      C = C + W(KK)
      PROFIT = PROFIT - P(KK)
      XX(KK) = XX(KK) - 1
      IF ( Z .LT. PROFIT + C*P(KK + 1)/W(KK + 1) ) GO TO 180
      C = C + XX(KK)*W(KK)
      PROFIT = PROFIT - XX(KK)*P(KK)
      XX(KK) = 0
      II = KK + 1
      GO TO 150
  180 IF ( R .LT. INT(MIN(KK)) ) GO TO 190
      II = KK + 1
      GO TO 60
  190 NN = KK + 1
      II = KK + 1
C
C STEP 7.
C
  200 DIFF = W(NN) - W(KK)
      IF ( DIFF ) 210,260,220
  210 T = R - DIFF
      IF ( T .LT. INT(MIN(NN-1)) ) GO TO 260
      S = C/W(NN)
      II = NN
      GO TO 70
  220 IF ( DIFF .GT. R ) GO TO 260
      IF ( Z .GE. PROFIT + P(NN) ) GO TO 260
      Z = PROFIT + P(NN)
      DO 230 J=1,KK
        X(J) = XX(J)
  230 CONTINUE
      KK1 = KK + 1
      DO 240 J=KK1,N
        X(J) = 0
  240 CONTINUE
      X(NN) = 1
      IF ( Z .NE. LIM ) GO TO 250
      C = CWF
      RETURN
  250 R = R - DIFF
      KK = NN
C
C STEP 8.
C
  260 IF ( Z .GE. PROFIT + C*P(NN+1)/W(NN+1) ) GO TO 150
      NN = NN + 1
      GO TO 200
      END
      SUBROUTINE REDU(N,PO,WO,JDIM,JPX,X)
C
C REDUCE AN UNBOUNDED KNAPSACK PROBLEM (PO,WO) THROUGH DOMINANCE
C RELATIONS.
C ON OUTPUT, JPX IS THE FIRST UNDOMINATED ITEM, X(JPX) THE SECOND, AND
C SO ON. IF Y IS THE LAST ONE, THEN X(Y) = 0.
C
      INTEGER PO(JDIM),WO(JDIM),X(JDIM)
      INTEGER FEQ,PRFEQ,PRJ,PRK
C INITIALIZATION.
      JPX = 1
      DO 10 J=1,N
        X(J) = J + 1
   10 CONTINUE
      X(N) = 0
      CRAT = FLOAT(PO(1))/FLOAT(WO(1)) + 1.
      PRFEQ = 0
      PRJ = 0
C MAIN ITERATION.
      J = JPX
   20   IWOJ = WO(J)
        IPOJ = PO(J)
        RJ = FLOAT(IPOJ)/FLOAT(IWOJ)
        IF ( RJ .EQ. CRAT ) GO TO 30
        CRAT = RJ
        PRFEQ = PRJ
        FEQ = J
        GO TO 80
C ITEMS K PRECEDING J WITH SAME RATIO.
   30   K = FEQ
        PRK = PRFEQ
   40     IF ( (WO(K)/IWOJ)*IPOJ .LT. PO(K) ) GO TO 60
C ITEM J DOMINATES ITEM K.
          IF ( PRK .EQ. 0 ) GO TO 50
          X(PRK) = X(K)
          GO TO 70
   50     JPX = X(K)
          GO TO 70
   60     PRK = K
   70     K = X(K)
        IF ( K .LT. J ) GO TO 40
C ITEMS K FOLLOWING J.
   80   K = X(J)
        PRK = J
   90     IF ( K .EQ. 0 ) GO TO 120
          IF ( (WO(K)/IWOJ)*IPOJ .LT. PO(K) ) GO TO 100
C ITEM J DOMINATES ITEM K.
          X(PRK) = X(K)
          GO TO 110
  100     PRK = K
  110     K = X(K)
        GO TO 90
  120   PRJ = J
        J = X(J)
      IF ( J .GT. 0 ) GO TO 20
      RETURN
      END
      SUBROUTINE SORTR(N,A,V,JDA)
C
C SORT THE REAL ARRAY A BY DECREASING VALUES (DERIVED FROM SUBROUTINE
C SORTZV OF THE C.E.R.N. LIBRARY).
C
C JDA           = LENGTH OF ARRAY A;
C N             = NUMBER OF ELEMENTS OF A TO BE SORTED;
C V(I) (INPUT)  = POINTER TO THE I-TH ELEMENT TO BE SORTED;
C V(I) (OUTPUT) = POINTER TO THE I-TH ELEMENT OF THE SORTED ARRAY.
C
C ON RETURN, ARRAY A IS UNCHANGED.
C
      INTEGER V(N),IU(20),IL(20)
      REAL    A(JDA)
      II = 1
      JJ = N
      IF ( N .LE. 1 ) RETURN
      M = 1
      I = II
      J = JJ
   10 IF ( I .GE. J ) GO TO 80
   20 K = I
      IJ = (J + I)/2
      IV = V(IJ)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 30
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
   30 L = J
      KI = V(J)
      IF ( A(KI) .LE. T ) GO TO 50
      V(IJ) = KI
      V(J) = IV
      IV = V(IJ)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 50
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
      GO TO 50
   40 V(L) = V(K)
      V(K) = IVT
   50 L = L - 1
      KI = V(L)
      IF ( A(KI) .LT. T ) GO TO 50
      IVT = KI
   60 K = K + 1
      KI = V(K)
      IF ( A(KI) .GT. T ) GO TO 60
      IF ( K .LE. L ) GO TO 40
      IF ( L - I .LE. J - K ) GO TO 70
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 90
   70 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 90
   80 M = M - 1
      IF ( M .EQ. 0 ) RETURN
      I = IL(M)
      J = IU(M)
   90 IF ( J - I .GE. II ) GO TO 20
      IF ( I .EQ. II ) GO TO 10
      I = I - 1
  100 I = I + 1
      IF ( I .EQ. J ) GO TO 80
      IV = V(I+1)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 100
      K = I
  110 V(K+1) = V(K)
      K = K - 1
      KI = V(K)
      IF ( T .GT. A(KI) ) GO TO 110
      V(K+1) = IV
      GO TO 100
      END
