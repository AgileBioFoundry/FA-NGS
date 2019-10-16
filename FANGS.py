import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

class Set_up:
    
    
### whole function needs the following sheets: source_plate, primers_destination, template_destination, and code.
    
    def __init__(self, file, plate_size):
        self.file = file
        self.plate_size =  plate_size 
    
        if plate_size == 96:
            columns = 12
            rows = 8
        elif plate_size == 384:
            columns = 24
            rows = 16
        converters = {col: str for col in range(columns + 1)}
        self.source_plate = pd.read_excel(file, 'source_plate', converters = converters)
        self.F_primer_dest = pd.read_excel(file, 'i7_dest_primers', converters = converters)
        self.R_primer_dest = pd.read_excel(file, 'i5_dest_primers', converters = converters)
        self.code = pd.read_excel(file, 'code')


    def Code(self):
        code = self.code
        
        loc = code[code.loc[:,"Transfer Volume"] == "Sample type"].index[0]    
        
        #resets code file for the rest of the function
        
        #tvol: transfer volume associated with a given sample type (e.g. primers are transfered at 150 nl)
        #part: defines the sample type of a Part Name (e.g. i501 is a primer)
        tvol_code = code[:loc].dropna(axis=0, how='any')
        part_code = code[loc+1:].dropna(axis=0, how='any')
        part_code.columns = ["Part Name","Sample type"]

        return  pd.merge(tvol_code, part_code)

    def alphanumerator(self, df, plate_type):
        #     label indeces
        
        df.index.name = 'rows'
        df.columns.name = 'columns'

        alphanumerated = df.unstack().dropna().reset_index()
        alphanumerated.rename(columns = {0:"Part Name"}, inplace =  True)
        alphanumerated = alphanumerated[alphanumerated.loc[:,"Part Name"] != " "]
        alphanumerated = alphanumerated.drop_duplicates(subset=['rows','columns']).reset_index(drop=True)

        
        if "source" in str(plate_type) :
            code =  self.Code()
            alphanumerated['Source Well'] = alphanumerated['rows'].map(str) + alphanumerated['columns'].map(str)
            alphanumerated = pd.merge(alphanumerated, code)

        elif "dest" in str(plate_type):
            alphanumerated['Destination Well'] = alphanumerated['rows'].map(str) + alphanumerated['columns'].map(str)

        else:
            print('''ERROR: plate_type can either be 'source' or 'destination', case sensitive''')


        return alphanumerated
        
    def PCR_set_up(self, num_columns = 24):
        
        set_up = Set_up(self.file, self.plate_size)



        indexed_F_primer_dest = set_up.alphanumerator(self.F_primer_dest, "dest")
        indexed_R_primer_dest = set_up.alphanumerator(self.R_primer_dest,"dest")
        indexed_source_plate = set_up.alphanumerator(self.source_plate, "source")



        F_prd = pd.merge(indexed_F_primer_dest,
                indexed_source_plate[['Part Name', 'Source Well', 'Transfer Volume', 'Sample type']], 
                         on='Part Name')
        R_prd = pd.merge(indexed_R_primer_dest,
                indexed_source_plate[['Part Name', 'Source Well', 'Transfer Volume', 'Sample type']],
                         on='Part Name')   

        echo_file = pd.concat([F_prd, R_prd])
        echo_file = echo_file.drop(['rows', 'columns'], axis=1)

        return echo_file





class Plates:
    
    def __init__(self, plate_size, quad= None, letters=None, numbers=None):
        self.plate_size = plate_size
        self.quad= quad


        
        if plate_size == 384:
            self.letters = 'ABCDEFGHIJKLMNOP'
            self.numbers = list(range(1,25))
        elif plate_size == 96:
            self.letters = 'ABCDEFGH'
            self.numbers = list(range(1,13))
        elif plate_size.lower() == "custom":
            self.letters = letters
            self.numbers = numbers

        else:
            print('''ERROR: plate_size can either be 96, 384 or  'custom' ''')
    
    def make_plate(self):
        plate = []
        num = self.numbers

        for letter in self.letters:
            if type(self.numbers) == int:
                num = list(range(1, self.numbers+1))
            for number in num:
                    plate.append('{0}{1}'.format(letter,number))

        return pd.DataFrame(plate)

    def make_quads(self):
        

        quad_map = Plates(384).make_plate()
        skip_rows = range(0,384,12)
         

        if self.quad == 1:
            quad_map = quad_map.iloc[0::2]
            quad_map = pd.concat([quad_map.iloc[skip_rows[0]:skip_rows[1]],
                                   quad_map.iloc[skip_rows[2]:skip_rows[3]],
                                   quad_map.iloc[skip_rows[4]:skip_rows[5]],
                                   quad_map.iloc[skip_rows[6]:skip_rows[7]],
                                   quad_map.iloc[skip_rows[8]:skip_rows[9]],
                                   quad_map.iloc[skip_rows[10]:skip_rows[11]],
                                   quad_map.iloc[skip_rows[12]:skip_rows[13]],
                                   quad_map.iloc[skip_rows[14]:skip_rows[15]]])            
        if self.quad == 2:
            quad_map = quad_map.iloc[1::2]
            quad_map = pd.concat([quad_map.iloc[skip_rows[0]:skip_rows[1]],
                                   quad_map.iloc[skip_rows[2]:skip_rows[3]],
                                   quad_map.iloc[skip_rows[4]:skip_rows[5]],
                                   quad_map.iloc[skip_rows[6]:skip_rows[7]],
                                   quad_map.iloc[skip_rows[8]:skip_rows[9]],
                                   quad_map.iloc[skip_rows[10]:skip_rows[11]],
                                   quad_map.iloc[skip_rows[12]:skip_rows[13]],
                                   quad_map.iloc[skip_rows[14]:skip_rows[15]]])   
        if self.quad == 3:
            quad_map = quad_map.iloc[0::2]
            quad_map = pd.concat([quad_map.iloc[skip_rows[1]:skip_rows[2]],
                                   quad_map.iloc[skip_rows[3]:skip_rows[4]],
                                   quad_map.iloc[skip_rows[5]:skip_rows[6]],
                                   quad_map.iloc[skip_rows[7]:skip_rows[8]],
                                   quad_map.iloc[skip_rows[9]:skip_rows[10]],
                                   quad_map.iloc[skip_rows[11]:skip_rows[12]],
                                   quad_map.iloc[skip_rows[13]:skip_rows[14]],
                                   quad_map.iloc[skip_rows[15]:skip_rows[16]]])   
        if self.quad == 4:
            quad_map = quad_map.iloc[1::2]
            quad_map = pd.concat([quad_map.iloc[skip_rows[1]:skip_rows[2]],
                                   quad_map.iloc[skip_rows[3]:skip_rows[4]],
                                   quad_map.iloc[skip_rows[5]:skip_rows[6]],
                                   quad_map.iloc[skip_rows[7]:skip_rows[8]],
                                   quad_map.iloc[skip_rows[9]:skip_rows[10]],
                                   quad_map.iloc[skip_rows[11]:skip_rows[12]],
                                   quad_map.iloc[skip_rows[13]:skip_rows[14]],
                                   quad_map.iloc[skip_rows[15]:skip_rows[16]]])

        quad_map = quad_map.rename(columns = {0:'quad_well'})
        quad_map = quad_map.reset_index().drop("index", axis=1)
        return quad_map
    
    
    
    def custom_plate_setup(self, failed_wells = None):

        if self.quad != None:
            plate = Plates(self.plate_size, self.quad, self.letters, self.numbers).make_quads().rename(columns = {"quad_well":"Source Well"})


        else:
            plate = Plates(self.plate_size, self.letters, self.numbers).make_plate().rename(columns = {0:"Source Well"})


        if failed_wells:
            print("failed wells:","\n",failed_wells )
            plate.set_index("Source Well", inplace=True)
            plate.drop(failed_wells, inplace = True)
            plate.reset_index(inplace = True)


        wells = list(plate.loc[:,"Source Well"])



        return wells
    
    def split_wells(self, failed_wells = None):
        wells = Plates(self.plate_size, self.quad, self.letters, self.numbers).custom_plate_setup(failed_wells= failed_wells)
        df = pd.DataFrame(wells).rename(columns = {0:"Wells"})
        well_ID = df.loc[:,'Wells'].str.split('(\d+)', expand=True).rename(columns = {0:"rows", 1:"columns"}).drop(2, axis=1)
        split_well = pd.concat([df, well_ID], axis =1 )
        
        
        return split_well




class Analysis:
    

    
    def __init__(self, file, plate_size, quad = None, letters = None, numbers = None, failed_wells = None):
        self.file = file
        self.plate_size = plate_size
        self.quad = quad
        self.letters = letters
        self.numbers = numbers
        self.failed_wells = failed_wells

        
    def read_deriv_file(self):
        
                    #read and edit derivative results file    

        df = pd.read_excel(self.file)     
        df.drop("Temperature", axis =1, inplace=True)
        df = df.transpose()
        
        wells = Plates(self.plate_size, self.quad, self.letters, self.numbers).custom_plate_setup(self.failed_wells)

        
        deriv_ranges = []
        for index in df.index.values:
            if index in wells:
                deriv_ranges.append(df.loc[index].values)
            else:
                pass
            
        return deriv_ranges
    
    def plot_melt_curve(self, failed_wells = None, color = "dodgerblue"):
        
        df = pd.read_excel(self.file)
        temp_range = np.array(df.loc[:,"Temperature"])
        plate = Plates(self.plate_size, self.quad, self.letters, self.numbers)
        deriv_ranges = Analysis(self.file, self.plate_size, self.quad, self.letters, self.numbers).read_deriv_file()
        wells = plate.custom_plate_setup()
        Well_ID = plate.split_wells()
        ncols = len(Well_ID.loc[:,"columns"].unique())
        nrows = len(Well_ID.loc[:,"rows"].unique())
        
        plt.style.use("default")


        fig,axes = plt.subplots(ncols= ncols,nrows= nrows, figsize = (2.5*ncols,2.5*nrows),
                                sharex=True,
                                sharey=True)
        plt.subplots_adjust(hspace=1)


        for ax, deriv_range, well in zip(axes.flat,deriv_ranges, wells):
            ax.plot(temp_range,deriv_range,color = color,linestyle = "-",linewidth=3)

            title = "{}".format(well)
            ax.set_title(title)

            ax = plt.gca() 
            
    def read_quant_file(self):
        df = pd.read_excel(self.file)     
        df.drop("Cycle", axis =1, inplace=True)
        df = df.transpose()
        
        wells = Plates(self.plate_size, self.quad, self.letters, self.numbers).custom_plate_setup(self.failed_wells)

        
        quant_ranges = []
        for index in df.index.values:
            if index in wells:
                quant_ranges.append(df.loc[index].values)
            else:
                pass
            
        return quant_ranges
    
    
    def plot_quant_curves(self, color = "dodgerblue"):
        
        df = pd.read_excel(self.file)     
        Cycles = np.array(df.loc[:,"Cycle"])
        plate = Plates(self.plate_size, self.letters, self.numbers)
        quant_ranges = Analysis(self.file, self.plate_size, self.quad, self.letters, self.numbers).read_quant_file()

        
        plt.style.use("default")


        fig,axes = plt.subplots(figsize = (3,2))
        plt.subplots_adjust(hspace=0.5)


        plt.hlines(1000, xmax=Cycles[-1], xmin = 0 ,
                   linestyle = "--", color = "grey", linewidth=.5 )

        for quant_range in quant_ranges:


            plt.plot(Cycles,quant_range,
                     color = color, linestyle = "-",
                     linewidth=1)
            ax = plt.gca()

        plt.xlabel("Cycle") 
        plt.ylabel("RFU")
        plt.ylim(top = 5000)
        plt.ylim(bottom = -100)
        
    def read_endRFU_file(self):
        split_plate = Plates(self.plate_size, self.quad, self.letters, self.numbers).split_wells()
        Well = []
        for row in split_plate.iterrows():
            if int(row[1].loc["columns"]) >=10:
                Well.append(row[1].rows+str(row[1].loc["columns"]))
            else:
                Well.append(row[1].rows+str(0)+str(row[1].loc["columns"]))
                
        split_plate["Well"] = Well
                
        df = pd.merge(pd.read_excel(self.file).reset_index().loc[:,['Well','End RFU']], split_plate).drop("Well", axis=1)        

        
        df_pivot = df.loc[:,['rows','columns','End RFU']].dropna(how='any').reset_index().pivot(index='rows',
                                                                                                columns='columns',
                                                                                                values='End RFU')
        return df_pivot
    
    def plot_end_RFU(self):
        df_pivot = Analysis(self.file, self.plate_size, self.quad, self.letters, self.numbers).read_endRFU_file()
        
        sns.set(context='paper', style='whitegrid', palette='deep', 
            font='sans-serif', font_scale=1.2, color_codes=True,  rc={"font.size":9})
        cols = len(df_pivot.columns.values)
        rows = len(df_pivot.index.values)

        if rows < 3:
            grid_kws = {"height_ratios": (.9, .4)}#, "hspace": 0.2*.7*rows}
            f, (ax, cbar_ax) = plt.subplots(2, figsize=(0.6*cols, 2*rows), gridspec_kw=grid_kws)

        else:
            grid_kws = {"height_ratios": (.9, .1)}#, "hspace": 0.2*.7*rows}
            f, (ax, cbar_ax) = plt.subplots(2, figsize=(0.6*cols, .4*rows), gridspec_kw=grid_kws)

        ax = sns.heatmap(df_pivot,  ax = ax, cbar_ax=cbar_ax,
                         linewidths = 0.05, 
                         cmap = sns.cubehelix_palette(n_colors=100, start=5.6, 
                                      rot=0.03, gamma=.3, hue=1.5, light=.8, dark=.01, reverse=False), 
                         vmin=0,
                         cbar_kws={'label': 'End_RFU', 
                                   "shrink": 0.2,
                                   "orientation": "horizontal"})

        ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 9)
        plt.tight_layout(pad=1)
        ax = plt.gca()


        
class Pooling_Calculator:
    def __init__(self, file, plate_size, quad = None, letters = None, numbers = None, failed_wells = None,
                max_volume = 500):
        self.file = file
        self.plate_size = plate_size
        self.quad = quad
        self.letters = letters
        self.numbers = numbers
        self.failed_wells = failed_wells
        self.max_volume = max_volume

    
    def read_end_RFU_file(self):
        
        split_plate = Plates(self.plate_size, self.quad, self.letters, self.numbers).split_wells(failed_wells = self.failed_wells)
        Well = []
        for row in split_plate.iterrows():
            if int(row[1].loc["columns"]) >=10:
                Well.append(row[1].rows+str(row[1].loc["columns"]))
            else:
                Well.append(row[1].rows+str(0)+str(row[1].loc["columns"]))
                
        split_plate["Well"] = Well
                
        df = pd.merge(pd.read_excel(self.file).reset_index().loc[:,['Well','End RFU']], split_plate).drop("Well", axis=1)        
        
        return df
    
    def make_pooling_file(self, max_transfer_vol = 15, dest_plate_size = 96):
    
        df = Pooling_Calculator(self.file, 
                                self.plate_size, 
                                self.quad, 
                                self.letters, 
                                self.numbers, 
                                self.failed_wells,
                                self.max_volume).read_end_RFU_file()
      
        min_End_RFU = df.loc[:,"End RFU"].min()
        
        df["Transfer Volume"] = self.max_volume*min_End_RFU/df.loc[:,"End RFU"]
        df.loc[:,"Transfer Volume"] =  df.loc[:,"Transfer Volume"].round(0)

            #Add the destination well(s) by minimizing total volume transfered per well

        df["Destination Well"] =""
        j = 0    
        i = 0
        transfer_wells = int(int(df.loc[:,"Transfer Volume"].sum()/1000)/max_transfer_vol)+1
        
        transfer_plate = Plates(plate_size= dest_plate_size).make_plate()

        for well in transfer_plate.iloc[:transfer_wells,0]:
            i = j

            if df.loc[i:,"Transfer Volume"].sum()/1000 < max_transfer_vol:
                df.loc[i:,"Destination Well"] = well

            else:  
                while df.loc[i:j,"Transfer Volume"].sum()/1000 < max_transfer_vol:
                    j = j + 1
                else:
                    df.loc[i:j,"Destination Well"] = well
        df.rename(columns = {"Wells":"Source Well"}, inplace=True)        
        return df
    
    
    def plot_transfer_volumes(self):
        
        df = df = Pooling_Calculator(self.file, 
                                self.plate_size, 
                                self.quad, 
                                self.letters, 
                                self.numbers, 
                                self.failed_wells,
                                self.max_volume).make_pooling_file()
        


        df_pivot = df.loc[:,['rows','columns','Transfer Volume']].dropna(how='any').reset_index().pivot(index='rows', 
                                                                               columns='columns', 
                                                                               values='Transfer Volume')
        
        
        
        #initialize plot
        sns.set(context='paper', style='whitegrid', palette='deep', 
                font='sans-serif', font_scale=1.2, color_codes=True,  rc={"font.size":9})
        cols = len(df_pivot.columns.values)
        rows = len(df_pivot.index.values)

        if rows < 3:
            grid_kws = {"height_ratios": (.9, .4)}#, "hspace": 0.2*.7*rows}
            f, (ax, cbar_ax) = plt.subplots(2, figsize=(0.6*cols, 2*rows), gridspec_kw=grid_kws)

        else:
            grid_kws = {"height_ratios": (.9, .1)}#, "hspace": 0.2*.7*rows}
            f, (ax, cbar_ax) = plt.subplots(2, figsize=(0.6*cols, .4*rows), gridspec_kw=grid_kws)

        ax = sns.heatmap(df_pivot,  ax = ax, cbar_ax=cbar_ax,
                         linewidths = 0.05, 
                         cmap = sns.cubehelix_palette(n_colors=100, start=5.6, 
                                      rot=0.03, gamma=.3, hue=1.5, light=.8, dark=.01, reverse=False), 
                         vmin=0,
                         cbar_kws={'label': 'Transfer Volume', 
                                   "shrink": 0.2,
                                   "orientation": "horizontal"})

        ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 9)
        plt.tight_layout(pad=1)
        ax = plt.gca()


