########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

import adodb,adodb_mxodbc,datetime

try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
        
class adodb_mssql(adodb_mxodbc.adodb_mxodbc):
    databaseType = 'mssql'
    dataProvider = 'mxodbc'
    sysDate = 'convert(datetime,convert(char,GetDate(),102),102)'
    sysTimeStamp = 'GetDate()'
    replaceQuote = "''"
    
    def _newcursor(self,rs):
        return cursor_mssql(rs,self)    
        
class cursor_mssql(adodb_mxodbc.cursor_mxodbc):
    def __init__(self,rs,conn):
        adodb_mxodbc.cursor_mxodbc.__init__(self,rs,conn)

    def Affected_Rows(self):
        return self._conn.GetOne('select @@rowcount')

    def Insert_ID(self):
        return self._conn.GetOne('select @@IDENTITY')

if __name__ == '__main__':
    db = adodb_mssql()
    db.Connect("PROVIDER=MSDASQL;DRIVER={SQL Server};SERVER=sherkhan;DATABASE=NorthWind;UID=adodb;PWD=natsoft;Trusted_Connection=No")
    adodb.Test(db)