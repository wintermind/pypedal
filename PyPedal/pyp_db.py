###############################################################################
# NAME: pyp_db.py
# VERSION: 2.0.2 (22JUNE2022)
# AUTHOR: John B. Cole, PhD (john.b.cole@gmail.com)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   connectToDatabase()
#   createDatabase()
#   createPedigreeTable()
#   deleteTable()
#   populatePedigreeTable()
#   doesTableExist()
#   tableCountRows()
#   tableDropRows()
###############################################################################
#   connectToDatabase_adodb()
#   createDatabase_adodb()
#   createPedigreeTable_adodb()
#   deleteTable_adodb()
#   populatePedigreeTable_adodb()
#   doesTableExist_adodb()
#   tableCountRows_adodb()
#   tableDropRows_adodb()
###############################################################################

## @package pyp_db
# pyp_db contains a set of procedures used to create, modify, and query pedigrees stored in relational databases.

import logging, math, os, string, sys
import pyp_io
import pyp_nrm
import pyp_utils

try:
    from PyPedal import adodb
except ImportError:
    print '[ERRROR]: Unable to import adodb in pyp_db.py!'

try:
    from PyPedal import sqlite3
except ImportError:
    print '[ERROR]: Unable to import sqlite3 in pyp_db.py!'
    logging.error('Unable to import sqlite3 in pyp_db.py!')

##
# connectToDatabase() opens a connection to a user-specified SQLite 3 database.
# @param pedobj A PyPedal pedigree object.
# @retval An sqlite3 connection on success, False otherwise.
def connectToDatabase(pedobj):
    """
    connectToDatabase() opens a connection to a user-specified SQLite 3 database
    (https://www.sqlite.org/index.html).
    """

    # Let's create the connection object.
    try:
        conn = sqlite3.connect(pedobj.kw['database_file'])
    except:
        print '[ERROR]: Unable to connect to SQLite 3 database %s in pyp_db/connectToDatabase()!' % \
              ( pedobj.kw['database_file'] )
        logging.error('Unable to connect to SQLite 3 database %s in pyp_db/connectToDatabase()!' %
                      pedobj.kw['database_file'] )
        conn = False

    # # conn = adodb.NewADOConnection(pedobj.kw['database_type'])
    # # SQLite
    # if pedobj.kw['database_type'] == 'sqlite':
    #     conn.Connect(database=pedobj.kw['database_name'])
    # # Postgres
    # elif pedobj.kw['database_type'] == 'postgres':
    #     conn.Connect(host=pedobj.kw['database_host'], \
    #         user=pedobj.kw['database_user'], \
    #         password=pedobj.kw['database_passwd'], \
    #         dbname=pedobj.kw['database_name'], \
    #         port=pedobj.kw['database_port'])
    # # MySQL
    # else:
    #     conn.Connect(pedobj.kw['database_host'], \
    #         pedobj.kw['database_user'], \
    #         pedobj.kw['database_passwd'], \
    #         pedobj.kw['database_name'])
    # if not conn:
    #     # If we can't connect to the specified database, try and create the database.
    #     # This will fail if the specified user does not have permission to create databases.
    #     created_table = createDatabaseTable(pedobj,conn)
    #     if created_table:
    #         conn = adodb.NewADOConnection(pedobj.kw['database_type'])
    #     else:
    #         conn = False

    # Return the connection (or False).
    return conn

##
# createPedigreeTable() creates a new pedigree table in a database.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @param drop Boolean indicating if the data should be dropped from an existing table with the same name
# @retval True on success, False otherwise.
def createPedigreeTable(pedobj, conn=False, drop=False):
    """
    createPedigreeTable() creates a new pedigree table in a database. Note that the table has a simple,
    fixed structure that may not include all attributes of a NewAnimal object.

    If the table already exists in the specified database a warning will be issued and the table will
    not be created. If the drop parameter is set to True then an existing table will be dropped and
    a new one created -- this may cause loss of data!
    """
    table_created = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/createPedigreeTable() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/createPedigreeTable() did not contain a valid connection and a connection to the database could not be established!')
    # If the table already exists we need to warn the user instead of clobbering their data.
    if doesTableExist(pedobj, conn=conn):
        if pedobj.kw['database_debug']:
            print '[INFO]: The table %s already exists in the database %s in pyp_db/createPedigreeTable().' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
        logging.info('The table %s already exists in the database %s in pyp_db/createPedigreeTable().', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Only delete data if specifically told to do so.
        if drop:
            tableDropRows(pedobj, conn=conn)
            if pedobj.kw['database_debug']:
                print '[WARNING]: Dropping rows from the table %s in the database %s in pyp_db/createPedigreeTable() because you told me to.' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.warning('Dropping rows from the table %s in the database %s in pyp_db/createPedigreeTable() because you told me to.', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Warn that data may be lost so the operation was cancelled.
        else:
            if pedobj.kw['database_debug']:
                print '[INFO]: The table %s in the database %s in pyp_db/createPedigreeTable() already exists and contains data that you did not tell me to delete.' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.info('The table %s in the database %s in pyp_db/createPedigreeTable() already exists and contains data that you did not tell me to delete.', pedobj.kw['database_table'], pedobj.kw['database_name'])
            if created_conn == True:
                conn.Close()
            return table_created
    # If you don't like the structure of the pedigree table then this is where you need to make changes. Make sure you
    # carefully check for side effects. For example, if you define a new table structure make sure you fix the loader
    # in pyp_newclasses so that your new attributes also are loaded.
    else:
        try:
            sql = 'CREATE TABLE %s (  \
                animalID INTEGER PRIMARY KEY, \
                animalName VARCHAR(128), \
                sireID INTEGER, \
                sireName VARCHAR(128), \
                damID INTEGER, \
                damName VARCHAR(128), \
                generation REAL, \
                infGeneration REAL, \
                birthyear INTEGER, \
                sex CHAR(1), \
                coi REAL, \
                founder CHAR(1), \
                ancestor CHAR(1), \
                originalID VARCHAR(128), \
                renumberedID INTEGER, \
                pedigreeComp REAL, \
                breed VARCHAR(128), \
                age REAL, \
                alive CHAR(1), \
                num_sons INTEGER, \
                num_daus INTEGER, \
                num_unk INTEGER, \
                herd INTEGER, \
                originalHerd VARCHAR(128), \
                gencoeff REAL, \
                alleles VARCHAR(256), \
                userField CHAR(128));' % ( pedobj.kw['database_table'] )
            #cursor = conn.Execute(sql)
            cursor = conn.cursor()
            cursor.execute(sql)
            conn.commit()
        # Crumbs...something went horribly wrong here!
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: Could not create the table %s in the database %s in pyp_db/connectToDatabase!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('Could not create the table %s in the database %s in pyp_db/connectToDatabase!', pedobj.kw['database_table'], pedobj.kw['database_name'])
    # Clean up any connections we opened. The database administrator will thank you for it.
    if created_conn == True:
        conn.close()
    return table_created

##
# deleteTable() drops a table from a database -- this can cause data loss if used carelessly!
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def deleteTable(pedobj, tablename=False, conn=False):
    """
    deleteTable() drops a table from a database, which can cause data loss if used carelessly!
    """
    # If no table name is specified use the default associated with the pedigree
    if tablename == False:
        tablename = pedobj.kw['database_table']
    table_dropped = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/deleteTable() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/deleteTable() did not contain a valid connection and a connection to the database could not be established!')
    # If we did get a good connection then let the massacre begin! Won't someone please think of the poor data?
    else:
        sql = 'DROP TABLE %s' % ( tablename )
        # Did it work? Yes!
        try:
            #cursor = conn.Execute(sql)
            cursor = conn.cursor()
            cursor.execute(sql)
            conn.commit()
            table_dropped = True
        # ...or not.
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: Could not delete the table %s from the database %s in pyp_db/deleteTable()!' % ( tablename, pedobj.kw['database_name'] )
            logging.error('Could not delete the table %s from the database %s in pyp_db/deleteTable()!', tablename, pedobj.kw['database_name'])
    # Clean up any connections we opened. The database administrator will thank you for it.
    if created_conn == True:
        conn.close()
    return table_dropped

##
# populatePedigreeTable() takes a PyPedal pedigree object and loads
# the animal records in that pedigree into a database table.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def populatePedigreeTable(pedobj,conn=False):
    """
    populatePedigreeTable() takes a PyPedal pedigree object and loads
    the animal records in that pedigree into a database table.
    """
    table_loaded = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/populatePedigreeTable() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/populatePedigreeTable() did not contain a valid connection and a connection to the database could not be established!')
    else:
        # If the pedigree table doesn't exist try and create it.
        if not doesTableExist(pedobj, conn=conn):
            #print 'Pedigree table does not exist!'
            created_table = createPedigreeTable(pedobj,conn)
            # If we can't create the table then we have to bail out.
            if not created_table:
                logging.error('Unable to create pedigree table in pyp_db/populatePedigreeTable()!')
            # Woohoo! We created the table!
            else:
                logging.info('Created pedigree table in pyp_db/populatePedigreeTable()!')
        # Okay, the pedigree table should now exist.
        if doesTableExist(pedobj, conn=conn):
            #print 'Pedigree table does exist!'
            try:
                for _p in pedobj.pedigree:
                    alleles = '__'.join(_p.alleles)
                    sql = 'INSERT INTO %s (animalID, animalName, sireID, sireName, damID,damName, generation, infGeneration, birthyear, sex, coi, founder, ancestor, originalID, renumberedID, pedigreeComp, breed, age, alive, num_sons, num_daus, num_unk, herd, gencoeff, alleles, userField) VALUES (%d, \'%s\', %d, \'%s\', %d, \'%s\', %s, %s, %d, \'%s\', %s, \'%s\', \'%s\', \'%s\', \'%s\', %s, \'%s\', %d, \'%s\', %d, %d, %d, \'%s\', %s, \'%s\', \'%s\' )' % ( pedobj.kw['database_table'], _p.animalID, _p.name, _p.sireID, _p.sireName, _p.damID, _p.damName, float(_p.gen), float(_p.igen), int(_p.by), _p.sex, float(_p.fa), _p.founder, _p.ancestor, _p.originalID, _p.originalID , float(_p.pedcomp), _p.breed, int(_p.age), _p.alive, len(_p.sons), len(_p.daus), len(_p.unks), _p.originalHerd, float(_p.gencoeff), alleles, str(_p.userField) )
                    # print sql
                    # cursor = conn.Execute(sql)
                    cursor = conn.cursor()
                    cursor.execute(sql)
                    conn.commit()
            except:
                if pedobj.kw['database_debug']:
                    print '[ERROR]: Unable to write to the table %s in the database %s!' % \
                        ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
                logging.error('Unable to write to the table %s in the database %s!', \
                    pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Try and create the table
        else:
            if pedobj.kw['database_debug']:
                print '[ERROR]: The pedigree table %s does not exist in the database %s in pyp_db/populatePedigreeTable()!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('The pedigree table %s does not exist in the database %s in pyp_db/populatePedigreeTable()!', pedobj.kw['database_table'], pedobj.kw['database_name'])
    if created_conn == True:
        conn.close()
        del conn
    return table_loaded

##
# doesTableExist() queries the database to determine whether or not a table exists.
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def doesTableExist(pedobj,tablename=False,conn=False):
    """
    doesTableExist() queries the database to determine whether or not a table exists.
    """
    if tablename == False:
        tablename = pedobj.kw['database_table']
    table_exists = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/doesTableExist() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/doesTableExist() did not contain a valid connection and a connection to the database could not be established!')
    else:
        try:
            sql = 'SELECT COUNT(*) FROM %s' % ( pedobj.kw['database_table'] )
            # cursor = conn.Execute(sql)
            cursor = conn.cursor()
            cursor.execute(sql)
            conn.commit()
            table_exists = True
            # cursor.Close()
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: The table %s does not exist in the database %s!' % \
                    ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('The table %s does not exist in the database %s!', \
                pedobj.kw['database_table'], pedobj.kw['database_name'])
    if created_conn == True:
        conn.close()
    return table_exists

##
# tableCountRows() returns the number of rows in a table.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @retval An integer on success, 0 otherwise
def tableCountRows(pedobj,conn=False):
    """
    tableCountRows() returns the number of rows in a table.
    """
    table_rows = 0
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/tableCountRows() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/tableCountRows() did not contain a valid connection and a connection to the database could not be established!')
    else:
        if doesTableExist(pedobj,conn=conn):
            try:
                sql = 'SELECT COUNT(animalID) FROM %s' % ( pedobj.kw['database_table'] )
                # cursor = conn.Execute(sql)
                cursor = conn.cursor()
                for row in cursor.execute(sql):
                    table_rows = row[0]
                # cursor.Close()
            except:
                table_rows = 0
        # doesTableExist() should have logged the table's non-existance for us, so we
        # can simply carry on.
        else:
            pass
    if created_conn == True:
        conn.close()
    return table_rows

##
# tableDropRows() deletes the rows from an existing table
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def tableDropRows(pedobj,tablename=False,conn=False):
    """
    tableDropRows() deletes the rows from an existing table
    """
    if tablename == False:
        tablename = pedobj.kw['database_table']
    rows_dropped = 0
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/tableDropRows() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/tableDropRows() did not contain a valid connection and a connection to the database could not be established!')
    else:
        if doesTableExist(pedobj,tablename,conn):
            try:
                sql = 'DELETE * FROM %s' % ( pedobj.kw['database_table'] )
                # cursor = conn.Execute(sql)
                cursor = conn.cursor()
                rows = cursor.execute(sql).rowcount
                conn.commit()
                rows_dropped = True
                # cursor.Close()
            except:
                if pedobj.kw['messages'] != 'quiet':
                    print '[ERROR]: pyp_db/tableDropRows() could not delete rows from the table %s in the database %s!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
                logging.error('pyp_db/tableDropRows() could not delete rows from the table %s in the database %s!', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # doesTableExist() should have logged the table's non-existance for us, so we
        # can simply carry on.
        else:
            pass
    if created_conn == True:
        conn.close()
    return rows_dropped

###############################################################################
# These are the old ADODB-based procedures for backwards compatibility.
###############################################################################

##
# connectToDatabase_adodb() opens a connection to a user-specified database.
# @param pedobj A PyPedal pedigree object.
# @retval An ADOdb connection on success, False otherwise.
def connectToDatabase_adodb(pedobj):
    """
    connectToDatabase_adodb() opens a connection to a user-specified database.
    """
    # These are the drivers recognized by SQLAlchemy's create_engine() URL Arguments
    # (see: http://www.sqlalchemy.org/docs/04/dbengine.html).
    drivers = ['mysql', 'postgres', 'sqlite']
    if pedobj.kw['database_type'] not in drivers:
        if pedobj.kw['database_debug']:
            print '[ERROR]: The database type %s is not recognized by pyp_db/connectToDatabase!' % \
                ( pedobj.kw['database_type'] )
        logging.error('The database type %s is not recognized by pyp_db/connectToDatabase!', \
            pedobj.kw['database_type'])
        conn = False

    # Let's create the connection object.
    conn = adodb.NewADOConnection(pedobj.kw['database_type'])
    # SQLite
    if pedobj.kw['database_type'] == 'sqlite':
        conn.Connect(database=pedobj.kw['database_name'])
    # Postgres
    elif pedobj.kw['database_type'] == 'postgres':
        conn.Connect(host=pedobj.kw['database_host'], \
            user=pedobj.kw['database_user'], \
            password=pedobj.kw['database_passwd'], \
            dbname=pedobj.kw['database_name'], \
            port=pedobj.kw['database_port'])
    # MySQL
    else:
        conn.Connect(pedobj.kw['database_host'], \
            pedobj.kw['database_user'], \
            pedobj.kw['database_passwd'], \
            pedobj.kw['database_name'])
    if not conn:
        # If we can't connect to the specified database, try and create the database.
        # This will fail if the specified user does not have permission to create databases.
        created_table = createDatabaseTable_adodb(pedobj,conn)
        if created_table:
            conn = adodb.NewADOConnection(pedobj.kw['database_type'])
        else:
            conn = False
    # Return the connection (or False).
    return conn

##
# createPedigreeDatabase_adodb() creates a new pedigree table in a database.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @param drop Boolean indicating if the data should be dropped from an existing table with the same name
# @retval True on success, False otherwise.
def createPedigreeTable_adodb(pedobj,conn=False,drop=False):
    """
    createPedigreeTable_adodb() creates a new pedigree table in a database. Note that the table has a simple,
    fixed structure that may not include all attributes of a NewAnimal object.

    If the table already exists in the specified database a warning will be issued and the table will
    not be created. If the drop parameter is set to True then an existing table will be dropped and
    a new one created -- this may cause loss of data!
    """
    table_created = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/createPedigreeTable_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/createPedigreeTable_adodb() did not contain a valid connection and a connection to the database could not be established!')
    # If the table already exists we need to warn the user instead of clobbering their data.
    if doesTableExist_adodb(pedobj,conn=conn):
        if pedobj.kw['database_debug']:
            print '[INFO]: The table %s already exists in the database %s in pyp_db/createPedigreeTable_adodb().' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
        logging.info('The table %s already exists in the database %s in pyp_db/createPedigreeTable_adodb().', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Only delete data if specifically told to do so.
        if drop:
            tableDropRows_adodb(pedobj,conn=conn)
            if pedobj.kw['database_debug']:
                print '[WARNING]: Dropping rows from the table %s in the database %s in pyp_db/createPedigreeTable_adodb() because you told me to.' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.warning('Dropping rows from the table %s in the database %s in pyp_db/createPedigreeTable_adodb() because you told me to.', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Warn that data may be lost so the operation was cancelled.
        else:
            if pedobj.kw['database_debug']:
                print '[INFO]: The table %s in the database %s in pyp_db/createPedigreeTable_adodb() already exists and contains data that you did not tell me to delete.' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.info('The table %s in the database %s in pyp_db/createPedigreeTable_adodb() already exists and contains data that you did not tell me to delete.', pedobj.kw['database_table'], pedobj.kw['database_name'])
            if created_conn == True:
                conn.Close()
            return table_created
    # If you don't like the structure of the pedigree table then this is where you need to make changes. Make sure you
    # carefully check for side effects. For example, if you define a new table structure make sure you fix the loader
    # in pyp_newclasses so that your new attributes also are loaded.
    else:
        try:
            sql = 'CREATE TABLE %s (  \
                animalID INTEGER PRIMARY KEY, \
                animalName VARCHAR(128), \
                sireID INTEGER, \
                sireName VARCHAR(128), \
                damID INTEGER, \
                damName VARCHAR(128), \
                generation REAL, \
                infGeneration REAL, \
                birthyear INTEGER, \
                sex CHAR(1), \
                coi REAL, \
                founder CHAR(1), \
                ancestor CHAR(1), \
                originalID VARCHAR(128), \
                renumberedID INTEGER, \
                pedigreeComp REAL, \
                breed VARCHAR(128), \
                age REAL, \
                alive CHAR(1), \
                num_sons INTEGER, \
                num_daus INTEGER, \
                num_unk INTEGER, \
                herd INTEGER, \
                originalHerd VARCHAR(128), \
                gencoeff REAL, \
                alleles VARCHAR(256), \
                userField CHAR(128));' % ( pedobj.kw['database_table'] )
            cursor = conn.Execute(sql)
            cursor.Close()
        # Crumbs...something went horribly wrong here!
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: Could not create the table %s in the database %s in pyp_db/connectToDatabase_adodb!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('Could not create the table %s in the database %s in pyp_db/connectToDatabase_adodb!', pedobj.kw['database_table'], pedobj.kw['database_name'])
    # Clean up any connections we opened. The database administrator will thank you for it.
    if created_conn == True:
        conn.Close()
    return table_created

##
# deleteTable_adodb() drops a table from a database -- this can cause data loss if used carelessly!
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def deleteTable_adodb(pedobj, tablename=False, conn=False):
    """
    deleteTable_adodb() drops a table from a database, which can cause data loss if used carelessly!
    """
    # If no table name is specified use the default associated with the pedigree
    if tablename == False:
        tablename = pedobj.kw['database_table']
    table_dropped = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/deleteTable_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/deleteTable_adodb() did not contain a valid connection and a connection to the database could not be established!')
    # If we did get a good connection then let the massacre begin! Won't someone please think of the poor data?
    else:
        sql = 'DROP TABLE %s' % ( tablename )
        # Did it work? Yes!
        try:
            cursor = conn.Execute(sql)
            table_dropped = True
            cursor.Close()
        # ...or not.
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: Could not delete the table %s from the database %s in pyp_db/deleteTable_adodb()!' % ( tablename, pedobj.kw['database_name'] )
            logging.error('Could not delete the table %s from the database %s in pyp_db/deleteTable_adodb()!', tablename, pedobj.kw['database_name'])
    # Clean up any connections we opened. The database administrator will thank you for it.
    if created_conn == True:
        conn.Close()
    return table_dropped

##
# populatePedigreeTable_adodb() takes a PyPedal pedigree object and loads
# the animal records in that pedigree into a database table.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def populatePedigreeTable_adodb(pedobj,conn=False):
    """
    populatePedigreeTable_adodb() takes a PyPedal pedigree object and loads
    the animal records in that pedigree into a database table.
    """
    table_loaded = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/populatePedigreeTable_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/populatePedigreeTable_adodb() did not contain a valid connection and a connection to the database could not be established!')
    else:
        # If the pedigree table doesn't exist try and create it.
        if not doesTableExist_adodb(pedobj, conn=conn):
            #print 'Pedigree table does not exist!'
            created_table = createPedigreeTable_adodb(pedobj,conn)
            # If we can't create the table then we have to bail out.
            if not created_table:
                logging.error('Unable to create pedigree table in pyp_db/populatePedigreeTable_adodb()!')
            # Woohoo! We created the table!
            else:
                logging.info('Created pedigree table in pyp_db/populatePedigreeTable_adodb()!')
        # Okay, the pedigree table should now exist.
        if doesTableExist_adodb(pedobj, conn=conn):
            #print 'Pedigree table does exist!'
            try:
                for _p in pedobj.pedigree:
                    alleles = '__'.join(_p.alleles)
                    sql = 'INSERT INTO %s (animalID, animalName, sireID, sireName, damID,damName, generation, infGeneration, birthyear, sex, coi, founder, ancestor, originalID, renumberedID, pedigreeComp, breed, age, alive, num_sons, num_daus, num_unk, herd, gencoeff, alleles, userField) VALUES (%d, \'%s\', %d, \'%s\', %d, \'%s\', %s, %s, %d, \'%s\', %s, \'%s\', \'%s\', \'%s\', \'%s\', %s, \'%s\', %d, \'%s\', %d, %d, %d, \'%s\', %s, \'%s\', \'%s\' )' % ( pedobj.kw['database_table'], _p.animalID, _p.name, _p.sireID, _p.sireName, _p.damID, _p.damName, float(_p.gen), float(_p.igen), int(_p.by), _p.sex, float(_p.fa), _p.founder, _p.ancestor, _p.originalID, _p.originalID , float(_p.pedcomp), _p.breed, int(_p.age), _p.alive, len(_p.sons), len(_p.daus), len(_p.unks), _p.originalHerd, float(_p.gencoeff), alleles, str(_p.userField) )
                    #print sql
                    cursor = conn.Execute(sql)
                    cursor.Close()
            except:
                if pedobj.kw['database_debug']:
                    print '[ERROR]: Unable to write to the table %s in the database %s!' % \
                        ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
                logging.error('Unable to write to the table %s in the database %s!', \
                    pedobj.kw['database_table'], pedobj.kw['database_name'])
        # Try and create the table
        else:
            if pedobj.kw['database_debug']:
                print '[ERROR]: The pedigree table %s does not exist in the database %s in pyp_db/populatePedigreeTable_adodb()!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('The pedigree table %s does not exist in the database %s in pyp_db/populatePedigreeTable_adodb()!', pedobj.kw['database_table'], pedobj.kw['database_name'])
    if created_conn == True:
        conn.Close()
        del conn
    return table_loaded

##
# doesTableExist_adodb() queries the database to determine whether or not a table exists.
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def doesTableExist_adodb(pedobj, tablename=False, conn=False):
    """
    doesTableExist() queries the database to determine whether or not a table exists.
    """
    if tablename == False:
        tablename = pedobj.kw['database_table']
    table_exists = False
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/doesTableExist_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/doesTableExist_adodb() did not contain a valid connection and a connection to the database could not be established!')
    else:
        try:
            sql = 'SELECT COUNT(*) FROM %s' % ( pedobj.kw['database_table'] )
            cursor = conn.Execute(sql)
            table_exists = True
            cursor.Close()
        except:
            if pedobj.kw['database_debug']:
                print '[ERROR]: The table %s does not exist in the database %s!' % \
                    ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
            logging.error('The table %s does not exist in the database %s!', \
                pedobj.kw['database_table'], pedobj.kw['database_name'])
    if created_conn == True:
        conn.Close()
    return table_exists

##
# tableCountRows_adodb() returns the number of rows in a table.
# @param pedobj A PyPedal pedigree object.
# @param conn An existing ADOdb connection or False to create one
# @retval An integer on success, 0 otherwise
def tableCountRows_adodb(pedobj, conn=False):
    """
    tableCountRows_adodb() returns the number of rows in a table.
    """
    table_rows = 0
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn = True
    else:
        logging.info('pyp_db/tableCountRows_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/tableCountRows_adodb() did not contain a valid connection and a connection to the database could not be established!')
    else:
        if doesTableExist(pedobj,conn=conn):
            try:
                sql = 'SELECT COUNT(animalID) FROM %s' % ( pedobj.kw['database_table'] )
                cursor = conn.Execute(sql)
                for row in cursor:
                    table_rows = row[0]
                cursor.Close()
            except:
                table_rows = 0
        # doesTableExist_adodb() should have logged the table's non-existance for us, so we
        # can simply carry on.
        else:
            pass
    if created_conn == True:
        conn.Close()
    return table_rows

##
# tableDropRows_adodb() deletes the rows from an existing table
# @param pedobj A PyPedal pedigree object.
# @param tablename The name of the table to delete.
# @param conn An existing ADOdb connection or False to create one
# @retval True on success, False otherwise.
def tableDropRows_adodb(pedobj,tablename=False,conn=False):
    """
    tableDropRows_adodb() deletes the rows from an existing table
    """
    if tablename == False:
        tablename = pedobj.kw['database_table']
    rows_dropped = 0
    created_conn = False
    # If the user doesn't pass us a conn then try and connect to the database
    if conn == False:
        conn = connectToDatabase_adodb(pedobj)
        created_conn == True
    else:
        logging.info('pyp_db/tableDropRows_adodb() established a connection to the database!')
    # If the user didn't give us a connection and we weren't able to connect to the
    # database ourselves then we have to give up.
    if conn == False:
        logging.error('The conn passed to pyp_db/tableDropRows_adodb() did not contain a valid connection and a connection to the database could not be established!')
    else:
        if doesTableExist_adodb(pedobj,tablename,conn):
            try:
                sql = 'DELETE * FROM %s' % ( pedobj.kw['database_table'] )
                cursor = conn.Execute(sql)
                rows = cursor.Affected_Rows()
                rows_dropped = True
                cursor.Close()
            except:
                if pedobj.kw['messages'] != 'quiet':
                    print '[ERROR]: pyp_db/tableDropRows_adodb() could not delete rows from the table %s in the database %s!' % ( pedobj.kw['database_table'], pedobj.kw['database_name'] )
                logging.error('pyp_db/tableDropRows_adodb() could not delete rows from the table %s in the database %s!', pedobj.kw['database_table'], pedobj.kw['database_name'])
        # doesTableExist_adodb() should have logged the table's non-existance for us, so we
        # can simply carry on.
        else:
            pass
    if created_conn == True:
        conn.Close()
    return rows_dropped