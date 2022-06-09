program owenst
    version 14
    syntax varlist
    plugin call COwensTFunc `varlist'
end


if "`c(os)'" == "Windows" {
    if "`c(osdtl)'" == "64-bit" {
	    program COwensTFunc, plugin using(owenst_plugin_win64.dll)
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
