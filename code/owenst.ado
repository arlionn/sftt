program owenst
    version 14
    syntax varlist
	plugin call COwensTFunc `varlist'
end


if "`c(os)'" == "Windows" {
    if "`c(osdtl)'" == "64-bit" {
	    program COwensTFunc, plugin using(COwensTFunc.plugin.win64)
	}
	else {
	    program COwensTFunc, plugin using(COwensTFunc.plugin.win32)
	}
}
else if "`c(os)'" == "MacOSX" {
	program COwensTFunc, plugin using(COwensTFunc.plugin.mac)
}
else {
    program COwensTFunc, plugin using(COwensTFunc.plugin.linux)  
}



// Windows: g++ -shared -DSYSTEM=STWIN32 -fPIC stplugin.c CSkewNormal.cpp -o COwensTFunc.plugin.win64 -I C:\boost_1_74_0
// macOS: g++ -bundle -DSYSTEM=APPLEMAC stplugin.c CSkewNormal.cpp -o COwensTFunc.plugin.mac
// Linux: g++ -shared -DSYSTEM=OPUNIX -fPIC stplugin.c CSkewNormal.cpp -o COwensTFunc.plugin.linux
