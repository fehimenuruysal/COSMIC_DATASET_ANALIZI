import pandas as pd

class TSV():

    def __init__(self):
        self.dosya_adi = "CosmicMutantExportCensus.tsv"
        self.dosya = pd.read_csv(self.dosya_adi , sep="\t")
        pd.set_option('display.max_rows' , self.dosya.shape[0]+1)

    def tsv_dosya_yaz(self,isim):
        yazilacak_dosya = open(r'C:\\Users\ASUS\Desktop\Uygulama_Odevi_2\\'+ isim +'.txt', "w")
        yazilacak_dosya.write(str(self.dosya))

    def analiz_dosya_yaz(self,sonuc,isim):
        yazilacak_dosya = open(r'C:\\Users\ASUS\Desktop\Uygulama_Odevi_2\\'+ isim +'.txt', "w")
        yazilacak_dosya.write(str(sonuc))

    def hasta_id_getir(self,primary_site , primary_histology,hastalik_adi):
        self.dosya = pd.read_csv(self.dosya_adi, sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        self.dosya.query("Primary_site == '" + primary_site + "'", inplace=True)
        self.dosya.query("Primary_histology == '" + primary_histology + "'", inplace=True)
        hastalik = self.dosya.drop_duplicates(subset=['ID_sample'])
        pd.set_option('display.max_rows', hastalik.shape[0] + 1)
        self.analiz_dosya_yaz(hastalik[['ID_sample']],hastalik_adi+"_hasta_id")
        return float(hastalik[['ID_sample']].count())

    def hasta_mutasyonlu_gen_listesi(self,primary_site,primary_histology,hastalik_adi):
        self.dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        self.dosya.query("Primary_site == '" + primary_site + "' & Primary_histology == '" + primary_histology + "'", inplace=True)
        mutasyonlu_gen_listesi = self.dosya[['ID_sample', 'GENOMIC_MUTATION_ID']]
        mutasyonlu_gen_listesi=mutasyonlu_gen_listesi.sort_values(by='ID_sample')
        self.analiz_dosya_yaz(mutasyonlu_gen_listesi,hastalik_adi+"_hastalarinin_mutasyonlu_genleri")

    def gen_bul(self,primary_site,primary_histology):
        self.dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        self.dosya.query("Primary_site == '" + primary_site + "'", inplace=True)
        self.dosya.query("Primary_histology == '" + primary_histology + "'", inplace=True)
        self.dosya = self.dosya.drop_duplicates(subset=['Gene_name'])
        sonuc = self.dosya[['Gene_name']]
        return sonuc

    def ortak_gen_bul(self,primary_site_1,primary_histology_1 ,primary_site_2,primary_histology_2):
        hastalik_1 = self.gen_bul(primary_site_1,primary_histology_1)
        hastalik_2 = self.gen_bul(primary_site_2, primary_histology_2)
        birlesim_liste = pd.concat([hastalik_1, hastalik_2])
        birlesim_liste["duplicated"] = birlesim_liste.duplicated()
        birlesim_liste.query("duplicated == True", inplace=True)
        self.analiz_dosya_yaz(birlesim_liste[['Gene_name']],"ortak_gen_" + primary_site_1 + primary_histology_2 + "&" + primary_site_2+primary_histology_2)

    def en_cok_gorunen_gen_yuzde(self,toplam_hasta_sayisi,primary_site,primary_histology ):
        dosya_1 = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        dosya_1.columns = [column.replace(" ", "_") for column in dosya_1.columns]
        dosya_1.query("Primary_site == '" + primary_site + "'", inplace=True)
        dosya_1.query("Primary_histology == '" + primary_histology + "'", inplace=True)
        dosya_1.drop_duplicates(subset=['ID_sample'])
        dosya_1 = dosya_1.groupby(['Gene_name'])['ID_sample'].count().reset_index(name="Hasta_Sayisi")

        dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        dosya.columns = [column.replace(" ", "_") for column in dosya.columns]
        dosya.query("Primary_site == '" + primary_site + "'", inplace=True)
        dosya.query("Primary_histology == '" + primary_histology + "'", inplace=True)
        dosya = dosya.groupby(['Gene_name'])['ID_sample'].count().sort_values(ascending=False).head(30).reset_index(name="Gorulme_Sayisi")

        dosya['Yüzde'] = dosya_1['Hasta_Sayisi'].mul(100)
        dosya['Yüzde'] = dosya['Yüzde'].divide(other=toplam_hasta_sayisi)
        self.analiz_dosya_yaz(dosya , primary_site + primary_histology + "_en_cok_gorunen_30_gen_ve_yuzdeleri")

    def hastalik_hasta_sayi_dagilimi(self):
        self.dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        self.dosya = self.dosya.drop_duplicates(subset=['ID_sample'])
        sonuc=self.dosya.groupby(['Primary_site'])['Primary_site'].count().sort_values(ascending=False).reset_index(name="Hasta_Sayisi")
        self.analiz_dosya_yaz(sonuc,"hastalik_hasta_sayisi")

    def hasta_mutasyon_sayisi(self):
        self.dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        sonuc = self.dosya.groupby(['ID_sample'])['ID_sample'].count().sort_values(ascending=False).reset_index(name="Mutasyon_Sayisi")
        self.analiz_dosya_yaz(sonuc,"hasta_mutasyon_sayilari")

    def hastalik_yas_ortalama(self):
        dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        dosya=dosya.drop_duplicates(subset=['ID_sample'])
        dosya=dosya.groupby(['Primary site'])['Age'].mean().reset_index(name="Yas_Ortalamasi")
        self.analiz_dosya_yaz(dosya,"hastalik_yas_ortalamalari")

    def mutasyonlu_gen_isimleri_gorulme_sayilari(self):
        dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        dosya = dosya.groupby(['GENOMIC_MUTATION_ID'])['ID_sample'].count().sort_values(ascending=False).reset_index(name="Görülme_Sayisi")
        self.analiz_dosya_yaz(dosya, "mutasyonlu_genler_ve_sayilari")

    def hastalik_yas_incelemesi(self,primary_site,primary_histology):
        self.dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        self.dosya.columns = [column.replace(" ", "_") for column in self.dosya.columns]
        self.dosya.query("Primary_site == '" + primary_site + "' & Primary_histology == '" + primary_histology + "'" ,inplace=True)
        sonuc ="En Yüksek Yaş : " + str(self.dosya['Age'].max())+"\nEn Düşük Yaş : " + str(self.dosya['Age'].min())+"\nOrtalama Yaş : " + str(self.dosya['Age'].mean())
        self.analiz_dosya_yaz(sonuc , primary_site+primary_histology+"_yas_incelemesi")

    def ortak_mutasyonlu_genler_3_hastalik(self):
        liver_dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        liver_dosya.columns = [column.replace(" ", "_") for column in liver_dosya.columns]
        liver_dosya.query("Primary_site == 'liver'", inplace=True)
        liver_dosya.query("Primary_histology == 'carcinoma'", inplace=True)
        liver = liver_dosya.drop_duplicates(subset=['GENOMIC_MUTATION_ID'])
        liver_listesi = liver[['GENOMIC_MUTATION_ID']]

        skin_dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        skin_dosya.columns = [column.replace(" ", "_") for column in skin_dosya.columns]
        skin_dosya.query("Primary_site == 'skin'", inplace=True)
        skin_dosya.query("Primary_histology == 'carcinoma'", inplace=True)
        skin = skin_dosya.drop_duplicates(subset=['GENOMIC_MUTATION_ID'])
        skin_listesi = skin[['GENOMIC_MUTATION_ID']]

        malignant_dosya = pd.read_csv("CosmicMutantExportCensus.tsv", sep="\t")
        malignant_dosya.columns = [column.replace(" ", "_") for column in malignant_dosya.columns]
        malignant_dosya.query("Primary_site == 'skin'", inplace=True)
        malignant_dosya.query("Primary_histology == 'malignant_melanoma'", inplace=True)
        malignant = malignant_dosya.drop_duplicates(subset=['GENOMIC_MUTATION_ID'])
        malignant_listesi = malignant[['GENOMIC_MUTATION_ID']]

        birlesim_liste = pd.concat([liver_listesi, skin_listesi])
        birlesim_liste["duplicated"] = birlesim_liste.duplicated()
        birlesim_liste.query("duplicated == True", inplace=True)
        birlesim_liste = birlesim_liste[['GENOMIC_MUTATION_ID']]

        birlesim_liste_2 = pd.concat([malignant_listesi, birlesim_liste])
        birlesim_liste_2["duplicated"] = birlesim_liste_2.duplicated()
        birlesim_liste_2.query("duplicated == True", inplace=True)
        pd.set_option('display.max_rows', birlesim_liste_2.shape[0] + 1)
        sonuc = birlesim_liste_2[['GENOMIC_MUTATION_ID']]
        self.analiz_dosya_yaz(sonuc , "malignant_skin_liver_ortak_mutasyonlu_genleri")

def main():
    analiz = TSV()
    analiz.tsv_dosya_yaz("tum_dosya")
    skin_sayisi = analiz.hasta_id_getir("skin", "carcinoma" , "skin_carcinoma")
    liver_sayisi = analiz.hasta_id_getir("liver", "carcinoma" , "liver_carcinom")
    analiz.hasta_mutasyonlu_gen_listesi("liver", "carcinoma" , "liver_carcinom")
    analiz.hasta_mutasyonlu_gen_listesi("skin", "carcinoma" , "skin_carcinoma")
    analiz.ortak_gen_bul("liver", "carcinoma" ,"skin", "carcinoma" )
    analiz.en_cok_gorunen_gen_yuzde(skin_sayisi,"skin", "carcinoma")
    analiz.en_cok_gorunen_gen_yuzde(liver_sayisi,"liver", "carcinoma")
    analiz.ortak_mutasyonlu_genler_3_hastalik()
    analiz.hastalik_hasta_sayi_dagilimi()
    analiz.hasta_mutasyon_sayisi()
    analiz.hastalik_yas_ortalama()
    analiz.mutasyonlu_gen_isimleri_gorulme_sayilari()
    analiz.hastalik_yas_incelemesi("liver", "carcinoma")
    analiz.hastalik_yas_incelemesi("skin", "carcinoma")

if __name__ == "__main__":
    main()