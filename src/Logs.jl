
"""������ ������� ���������� �� ������� � � ���-����. ������������, ��� logFileName � path ���������� �� ������������� �������"""
macro showAndLog(str)
    logfile = open(string(path, logFileName), "a+")
    str = string(eval(str), "\n")
    write(logfile, str)
    showall(str)
    close(logfile)
end
