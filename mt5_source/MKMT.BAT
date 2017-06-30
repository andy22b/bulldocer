@if %QB4%.==. set qb4=C:\BASIC\QB4
@if not exist %QB4%\qb.lib goto ndpth
@if %1.==0. goto dosub
@if exist m5sub.obj goto chkeig
:dosub
@set ON=m5sub
@if exist m5sub.obj del m5sub.obj
%QB4%\bc m5sub,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
goto doeig
:chkeig
@if exist m5eigen.obj goto chkpar
:doeig
@set ON=m5eigen
@if exist m5eigen.obj del m5eigen.obj
%QB4%\bc m5eigen,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
@if %2.==0. goto end
:chkpar
@if %1.==. goto 1
@if %1.==a. goto 7
@if %1.==A. goto 7
@echo off
FOR %%P IN (0 1 2 3 4 5 6 7 8 9 10 11) DO IF %1.==%%P. set ON=%1
echo on
@IF NOT %ON%.==. GOTO %ON%
@goto syntax
:0
:1
@set on=mt5men1
@if exist mt5men1.obj goto l1
%QB4%\bc mt5men1,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l1
link /e mt5men1+m5sub+m5eigen,,nul,%QB4%\qb+%QB4%\
@del mt5men1.obj
@if errorlevel 1 goto err
@if %2.==1. goto end
:2
@set on=mt5men2
@if exist mt5men2.obj goto l2
%QB4%\bc mt5men2,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l2
link /e mt5men2+m5sub+m5eigen,,nul,%QB4%\qb+%QB4%\
@del mt5men2.obj
@if errorlevel 1 goto err
@if %2.==2. goto end
:3
@set ON=m5sub3
@if exist m5sub3.obj goto 33
%QB4%\bc m5sub3,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:33
@set on=mt5it 
@if exist mt5it.obj goto l3
%QB4%\bc mt5it,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l3
link /e mt5it+m5sub+m5eigen+m5sub3,,nul,%QB4%\qb+%QB4%\
@del mt5it.obj
@if errorlevel 1 goto err
@if %2.==3. goto end
:4
@set ON=m5mtsub
@if exist m5mtsub.obj goto 45
%QB4%\bc m5mtsub,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:45
@set on=mt5plot
@if exist mt5plot.obj goto l4
%QB4%\bc mt5plot,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l4
link /e mt5plot+m5sub+m5mtsub,,nul,%QB4%\qb+%QB4%\
@del mt5plot.obj
@if errorlevel 1 goto err
@if %2.==4. goto end
:5
@set ON=m5sub1
@if exist m5sub1.obj goto 55
%QB4%\bc m5sub1,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:55
@set on=mt5main
@if exist mt5main.obj goto l5
%QB4%\bc mt5main,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l5
link /e mt5main+m5sub+m5sub1,,nul,%QB4%\qb+%QB4%\
@del mt5main.obj
@if errorlevel 1 goto err
@if %2.==5. goto end
:6
@set on=mt5
@if exist mt5.obj goto l6
%QB4%\bc mt5,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l6
link /e mt5+m5sub,,nul,%QB4%\qb+%QB4%\
@del mt5.obj
@if errorlevel 1 goto err
@if %2.==6. goto end
@echo Successful compilation of mt5 main programs!
@echo Hit Ctrl-Brk to ABORT compilation of auxiliary programs, or..
@set on=
@pause
:7
@set ON=m5sub7
@if exist m5sub7.obj goto 77
%QB4%\bc m5sub7,,nul /s /ah
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:77
@set on=mt5int
@if exist mt5int.obj goto l7
%QB4%\bc mt5int,,nul /s 
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l7
link /e mt5int+m5sub+m5sub7,,nul,%QB4%\qb+%QB4%\
@del mt5int.obj
@if errorlevel 1 goto err
@if %2.==7. goto end
:8
@set on=mt5gdsn
@if exist mt5gdsn.obj goto l8
%QB4%\bc mt5gdsn,,nul /s 
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l8
link /e mt5gdsn+m5sub,,nul,%QB4%\qb+%QB4%\
@del mt5gdsn.obj
@if errorlevel 1 goto err
@if %2.==8. goto end
:9
@set on=mt5list
@if exist mt5list.obj goto l9
%QB4%\bc mt5list,,nul /s 
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l9
link /e mt5list+m5sub,,nul,%QB4%\qb+%QB4%\
@del mt5list.obj
@if errorlevel 1 goto err
@if %2.==9. goto end
:10
@set on=mt5dig
@if exist mt5dig.obj goto l10
%QB4%\bc mt5dig,,nul /s
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l10
link /e mt5dig+m5sub,,nul,%QB4%\qb+%QB4%\
@del mt5dig.obj
@if errorlevel 1 goto err
@if %2.==10. goto end
:11
@set on=mt5prev
@if exist mt5prev.obj goto l11
%QB4%\bc mt5prev,,nul /s
@if errorlevel 1 del %on%.obj
@if errorlevel 1 goto err
:l11
link /e mt5prev+m5sub,,nul,%QB4%\qb+%QB4%\
@del mt5prev.obj
@if errorlevel 1 goto err
@if %2.==11. goto end
@echo Successful compilation of mt5 auxiliary programs!
@goto fin
:err
@echo Error occurred in module %ON%. Compilation terminated early.
@goto show
:ndpth
@echo Cannot locate %QB4%\QB.LIB in path %QB4%. Please issue command
@echo  SET QB4={....}  
@echo      where: {....} is full path to QB4 files (without '\' on end)
@echo             (default path is C:\BASIC\QB4, use . for current path)
@echo  Then run MKMT again.
@goto end
:syntax
@echo Syntax:
@echo  MKMT (para1 (para2))
@echo        where: () means optional (no entry means ALL)
@echo               para1 is starting module # (0-9,A=auxiliary mt5 files only)
@echo               para2 is ending module # (0-9, no entry means ALL REST)
:show
@echo  module #'s:
@echo   m5sub=0, mt5men1=1, mt5men2=2, mt5it=3, mt5plot=4, mt5main=5, mt5=6
@echo   Auxiliary modules: mt5int=7, mt5gdsn=8, mt5list=9, mt5dig=10, mt5prev=11
@goto fin
:end
@echo Compilation successful!
:fin
@set ON=
