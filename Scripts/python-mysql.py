  ## Connecting to the database

  ## importing 'mysql.connector' as mysql for convenient
  import mysql.connector as mysql

  ## connecting to the database using 'connect()' method
  ## it takes 3 required parameters 'host', 'user', 'passwd'
  db = mysql.connect(
      host = "localhost",
      user = "root",
      passwd = "lucas888"
  )

  print(db) # it will print a connection object if everything is fine


  ## We need to create a "cursor" object to navigate MySQL
  cursor = db.cursor()

  ## Create a DataBase
  cursor.execute("CREATE DATABASE PhD_Project")

  ## executing the statement using 'execute()' method
  cursor.execute("SHOW DATABASES")

  ## 'fetchall()' method fetches all the rows from the last executed statement
  databases = cursor.fetchall() ## it returns a list of all databases present

  ## printing the list of databases
  print(databases)

  cursor.close()
  db.close()

  filePath = "/media/lucas/Disc4T/Projects/PhD_Project/R_results/Stz_Areas/areaStz_anastEpiReset1.csv"

  with open(filePath, "r+") as f:
      header = f.readline().strip().replace("\"","").split(",")

  cols = [h.split("-")[0] + " FLOAT NULL" for h in header[:-1]]

  cmd = ", ".join(cols)
  cmd = cmd + ", Gene_Id VARCHAR(255) NOT NULL PRIMARY KEY"
  cmd = "CREATE TABLE IF NOT EXISTS aEpi1 ({});" .format(cmd)

  cmd

  db = mysql.connect(
        host = "localhost",
        user = "root",
        passwd = "lucas888",
        database="PhD_Project"
  )

  cursor = db.cursor()

  cursor.execute(cmd)
  #cursor.execute("DROP TABLE aEpi1")

  cursor.execute("SHOW TABLES")
  cursor.fetchall()

  cmd =

  CREATE TABLE users (id INT NOT NULL AUTO_INCREMENT, first_name VARCHAR(255) NOT NULL, last_name VARCHAR(255) NOT NULL, email VARCHAR(255) NOT NULL, transactions INT NOT NULL, account_creation DATE NOT NULL, PRIMARY KEY (id));

  LOAD DATA LOCAL INFILE '/Users/miguelgomez/Desktop/mock_data.csv' INTO TABLE users FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, first_name, last_name, email, transactions, @account_creation)SET account_creation  = STR_TO_DATE(@account_creation, '%m/%d/%y')

  def dict_from_CSV_lists(csvfile):

      """Read a csv file (GSEA format) and convert
      it to a dictionary with each header as a key
      and the rest of the column (as a list) as value."""

      with open(csvfile, "r+") as csv:
          first = True #On first iteration get dictionary keys
          for line in csv:
              if first:
                  keys = line.strip().split(",")
                  gene_lists_dict = {k:[] for k in keys}
                  first = False
              else:
                  genes = line.strip().split(",")
                  for i, key in enumerate(keys):
                      if genes[i] and genes[i] != "NA": # Check it is not an empty string or NA.
                          gene_lists_dict[key].append(genes[i])

      return(gene_lists_dict)


  def dict_from_GFF(gffile):

      """Read a GFF file and store it's info into a dict"""

      with open(gffile, "r+") as gff:
          gffdict = {}
          for line in gff:
              linelist = line.strip().split("\t")

              chrom = linelist[0]
              start = linelist[3]
              stop = linelist[4]
              pos = chrom+":"+start+"-"+stop

              strand = linelist[6]

              info = linelist[8]
              gid = info.split(";")[0].replace("ID=", "")
              annot = info.split(";")[1].replace("description=", "")

              gffdict[gid] = {"chrom":chrom,
                              "start":start,
                              "stop":stop,
                              "strand":strand,
                              "annot":annot}

      return(gffdict)


  csvfile = "/media/lucas/Disc4T/Projects/PhD_Project/External_Data/all_gene_lists_160719.csv"
  listdict = dict_from_CSV_lists(csvfile)

  ## There are repeated genes in the lists, let's remove them.
  for key, val in listdict.items():
      listdict[key] = list(set(val))


  gffile = "/media/lucas/Disc4T/Projects/Miniprojects/gene_gff.gff"
  gffdict = dict_from_GFF(gffile)

  import mysql.connector as mysql

  db = mysql.connect(
      host = "localhost",
      user = "root",
      passwd = "lucas888",
      database="PhD_Project"
  )

  print(db) # it will print a connection object if everything is fine
  cursor = db.cursor()

  cursor.execute("DROP TABLE gene_to_list")
  cursor.execute("DROP TABLE genes")
  cursor.execute("DROP TABLE geneLists")

  cmd = ("CREATE TABLE genes("
         "ID INT AUTO_INCREMENT PRIMARY KEY, "
         "plasmoID VARCHAR(255) UNIQUE, "
         "chrom VARCHAR(255), "
         "start INT, "
         "stop INT, "
         "strand CHAR(1), "
         "annot TEXT"
         ");")

  cursor.execute(cmd)

  cmd = ("CREATE TABLE geneLists("
         "ID INT AUTO_INCREMENT PRIMARY KEY, "
         "name VARCHAR(255) UNIQUE, "
         "source TEXT, "
         "description TEXT"
         ");")

  cursor.execute(cmd)

  cmd = ("CREATE TABLE gene_to_list("
         "geneID INT NOT NULL, "
         "listID INT NOT NULL, "
         "gene_plasmoID VARCHAR(255), "
         "list_name VARCHAR(255), "
         "FOREIGN KEY (geneID) REFERENCES genes (ID) "
         "ON DELETE RESTRICT ON UPDATE CASCADE, "
         "FOREIGN KEY (listID) REFERENCES geneLists (ID) "
         "ON DELETE RESTRICT ON UPDATE CASCADE, "
         "PRIMARY KEY (geneID, listID)"
         ");")

  cursor.execute(cmd)
  cursor.execute("SHOW TABLES")
  cursor.fetchall()
  cursor.close()
  db.close()

  import mysql.connector as mysql
  from mysql.connector import Error

  db = mysql.connect(
        host = "localhost",
        user = "root",
        passwd = "lucas888",
        database="PhD_Project"
  )
  cursor = db.cursor()

  ## Populate genes using gffdict

  cmd = ("INSERT INTO genes (plasmoID, chrom, start, stop, strand, annot) "
         "VALUES ('{}', '{d[chrom]}', {d[start]}, {d[stop]}, '{d[strand]}', \"{d[annot]}\");")

  for key, val in gffdict.items():
      try:
          insert = cmd .format(key, d = val)
          cursor.execute(insert)
      except mysql.connector.Error as error:
          print("Failed to insert {}: {}" .format(key, error))

  db.commit()

  ## Populate geneLists with some examples

  add_geneList = ("INSERT INTO geneLists (name, source, description) "
                  "VALUES ('{}', \"{}\", \"{}\");")

  cursor.execute(add_geneList .format("Fraschka_High_Confidence",
                                      "PMID:29649445",
                                      "List of HP1 enriched genes by ChIP-Seq."))

  cursor.execute(add_geneList .format("Lopez-Barragan",
                                      "PMID:22129310",
                                      "Mature gametocyte markers."))

  cursor.execute(add_geneList .format("Lasonder",
                                      "PMID:27298255",
                                      "Mature gametocyte markers."))

  db.commit()
  ## Create the rest of lists with same name as the respective
  ## column header and empty source and description.

  for k in listdict.keys():
      try:
          cursor.execute(add_geneList .format(k, "", ""))
          print("Created {}" .format(k))
      except mysql.Error as error:
          print("Failed to create {}: {}" .format(k, error))

  db.commit()


  ## Populate gene_to_list

  for key, value in listdict.items():

      get_geneID = "SELECT ID FROM genes WHERE plasmoID = '{}'"
      get_listID = "SELECT ID FROM geneLists WHERE name = '{}'"

      add_gene_to_list = ("INSERT INTO gene_to_list (geneID, listID, gene_plasmoID, list_name) "
                          "VALUES ({}, {}, '{}', '{}')")

      cursor.execute(get_listID .format(key))
      listID = cursor.fetchall()[0][0]

      for gen in value:

          cursor.execute(get_geneID .format(gen))
          answ = cursor.fetchall()
          try:
              geneID = answ[0][0]

          except:
              print(gen, answ)

          else:
              try:
                  cursor.execute(add_gene_to_list .format(geneID, listID, gen, key))
              except mysql.connector.Error as error:
                  print("Failed to insert {} {}: {}" .format(key, gen, error))

  db.commit()
  cursor.close()
  db.close()

  import mysql.connector as mysql
  from mysql.connector import Error

  db = mysql.connect(
      host = "localhost",
      user = "root",
      passwd = "lucas888",
      database="PhD_Project"
  )
  cursor = db.cursor()

  cursor.execute("DESCRIBE genes")
  cursor.fetchall()
  cursor.execute("DESCRIBE geneLists")
  cursor.fetchall()
  cursor.execute("DESCRIBE gene_to_list")
  cursor.fetchall()

  select = ("SELECT {} FROM {} WHERE {};")

  ## Query 1: Gene data from Genes in list

  query1 = ["g.plasmoID, g.annot, g.chrom, gtl.list_name",
            "genes as g INNER JOIN gene_to_list as gtl ON g.plasmoID = gtl.gene_plasmoID",
            "gtl.list_name = 'var'"]

  cursor.execute(select .format(*query1))
  hits = cursor.fetchall()

  len(hits)

  genes = []

  for pid, anot, chrom, lname in hits:
      print("{}\t{}\t{}\t{}" .format(pid, anot, chrom, lname))
      genes.append(pid)

  ## Query 2: Number of genes in every list

  select = ("SELECT {} FROM {};")


  query2 = ["COUNT(gene_plasmoID), list_name",
            "gene_to_list GROUP BY list_name"]

  print(select .format(*query2))


  cursor.execute(select .format(*query2))
  hits = cursor.fetchall()

  for h in hits:
      print("{}, {}" .format(h[1], h[0]))

  cursor.close()
  db.close()

  import pandas as pd

  ## Import ChIP-Seq data using pandas
  ## Import gene-wise coverage

  geneChIP = pd.read_csv(
      "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/GenWise_Coverage/all_coverage_RPKM.bed", sep="\t")

  geneChIP.head()

  geneChIP.columns = ['chrom', 'start', 'stop',
                      'gene_id', 'E5K9_me',
                      'NF54_ac', 'B11_in',
                      'C2_me', 'E5HA_in',
                      '1.2B_me', 'C2_in',
                      '3D7_in', 'B11_me',
                      '10G_me', '1.2B_in',
                      'B11_ac', 'NF54_in',
                      '10G_ac', 'A7K9_in',
                      'C2_ac', 'E5K9_in',
                      '1.2B_ac', 'E5HA_ac',
                      '10G_in', 'A7K9_me',
                      '3D7_me', 'NF54_me',
                      '3D7_ac', 'E5HA_me']

  promChIP = pd.read_csv(
     "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/GenWise_Coverage/all_coveragePromoters_RPKM.bed", sep="\t")

  promChIP = promChIP.drop(['chrom', 'start', 'stop'], axis=1)

  promChIP["gene_id"] = [g.replace("_promoter", "") for g in promChIP["gene_id"]]

  promChIP.columns = ['gene_id', 'E5K9_me_prom',
                      'NF54_ac_prom', 'B11_in_prom',
                      'C2_me_prom', 'E5HA_in_prom',
                      '1.2B_me_prom', 'C2_in_prom',
                      '3D7_in_prom', 'B11_me_prom',
                      '10G_me_prom', '1.2B_in_prom',
                      'B11_ac_prom', 'NF54_in_prom',
                      '10G_ac_prom', 'A7K9_in_prom',
                      'C2_ac_prom', 'E5K9_in_prom',
                      '1.2B_ac_prom', 'E5HA_ac_prom',
                      '10G_in_prom', 'A7K9_me_prom',
                      '3D7_me_prom', 'NF54_me_prom',
                      '3D7_ac_prom', 'E5HA_me_prom']

  allChip = geneChIP.join(promChIP.set_index('gene_id'), on='gene_id')

  # Set gene_id column as index and transpose df
  # (.to_dict) sets column names as keys and we want
  # gene_ids to be the keys.

  ## Standarize (-mean/std) each column
  ## Pandas applies operations on dfs columnwise by default.

  chipVals = allChip.iloc[:,4:]
  chipVals_stz = (chipVals-chipVals.mean())/chipVals.std()

  allChip.iloc[:,4:] = chipVals_stz
  allChip.head()

  chipDict = allChip.set_index("gene_id").T.to_dict("list")

  import mysql.connector as mysql
  from mysql.connector import Error

  db = mysql.connect(
      host="localhost",
      user="root",
      passwd="lucas888",
      database="PhD_Project"
  )
  cursor = db.cursor()

  ## Create chipSeq table

  #cursor.execute("DROP TABLE chipSeq")

  cmd = ("CREATE TABLE chipSeq("
         "ID INT AUTO_INCREMENT PRIMARY KEY, "
         "plasmoID VARCHAR(255), "
         "chrom VARCHAR(255), "
         "start INT NOT NULL, "
         "stop INT NOT NULL, "
         "{}"
         "FOREIGN KEY (ID) REFERENCES genes (ID) "
         "ON DELETE RESTRICT ON UPDATE CASCADE"
         ");")

  colstring = ""
  for col in allChip.columns[4:]:
      colname = col.replace(".","") #can not put "." in a column name (1.2B).
      colstring += "{} FLOAT, " .format(colname)

  cursor.execute(cmd .format(colstring))
  cursor.execute("SHOW TABLES")
  cursor.fetchall()

  # Populate table.

  cols = allChip.columns[4:]
  cols = [x.replace(".", "") for x in cols]
  cols = ", ".join(cols)

  cmd = ("INSERT INTO chipSeq (plasmoID, chrom, start, stop, {}) "
         "VALUES ('{}', '{}', {});")

  for key, val in chipDict.items():
      fields = ", ".join([str(v) for v in val[1:]])
      chrom = val[0]
      insert = cmd .format(cols, key, chrom, fields)
      #print(insert)
      try:
          cursor.execute(insert)

      except mysql.Error as error:
            print("Failed to insert {}: {}" .format(key, error))
            print("-------------")

  db.commit()

  cursor.close()
  db.close()

  import pandas as pd
  import mysql.connector as mysql
  from mysql.connector import Error

  db = mysql.connect(
      host="localhost",
      user="root",
      passwd="lucas888",
      database="PhD_Project"
  )
  cursor = db.cursor()

  ## Create chipSeq table

  cursor.execute("DROP TABLE variantome")

  cmd = ("CREATE TABLE variantome("
         "ID INT AUTO_INCREMENT PRIMARY KEY, "
         "plasmoID VARCHAR(255), "
         "10G FLOAT, "
         "12B FLOAT, "
         "E5K9 FLOAT, "
         "A7K9 FLOAT, "
         "B11 FLOAT, "
         "FOREIGN KEY (ID) REFERENCES genes (ID) "
         "ON DELETE RESTRICT ON UPDATE CASCADE"
         ");")

  cursor.execute(cmd)

  cursor.execute("SHOW TABLES")
  cursor.fetchall()

  # Populate table.

  trans = pd.read_csv("./External_Data/3D7variantome_10g12b3d7b_trans.csv", usecols=[1,2,3,4])
  transNULL

  trans["10G"] = pd.to_numeric(trans["10G"], errors="coerce")
  trans["1.2B"] = pd.to_numeric(trans["1.2B"], errors="coerce")
  trans["3D7-B"] = pd.to_numeric(trans["3D7-B"], errors="coerce")
  transNULL = trans.where((pd.notnull(trans)), "NULL")


  cmd = ("INSERT INTO variantome (plasmoID, 10G, 12B, E5K9, A7K9, B11) "
         "VALUES ('{}', {}, {}, {}, {}, {});")

  for index, row in transNULL.iterrows():
      insert = cmd .format(row[3], row[1], row[0], row[2], row[2], row[2])
      #print(insert)
      try:
          cursor.execute(insert)

      except mysql.Error as error:
          print("Failed to insert {}: {}" .format(row[3], error))
          print("-------------")

  db.commit()

  cursor.close()
  db.close()
