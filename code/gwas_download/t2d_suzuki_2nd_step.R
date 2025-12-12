# ---------------------------------------------------------------------------------- #
# 0 - K??t??phaneler
# ---------------------------------------------------------------------------------- #
library(tidyverse)
library(data.table)
library(readxl)
library(biomaRt) 

# (Burada senin dosya indirme kodlar??n zaten ??al??????yor, o k??sm?? ge??iyorum)
# Varsayal??m ki t2d_suzuki_lipo_st6 ve t2d_suzuki_gwas ??u an workspace'te y??kl??.

# ---------------------------------------------------------------------------------- #
# 1 - LIPO DATA STANDARDIZASYONU (t2d_suzuki_lipo_st6)
# ---------------------------------------------------------------------------------- #
# Kolonlar??n: chromosome, base_pair_location, variant, cluster_assignment...

lipo_clean <- t2d_suzuki_lipo_st6 %>%
  mutate(
    # String'e ??eviriyoruz ki say??/yaz?? kar??????kl?????? olmas??n
    chromosome = as.character(chromosome),
    base_pair_location = as.character(base_pair_location),
    # "chr" tak??s?? varsa siliyoruz (??rn: "chr1" -> "1")
    chromosome = gsub("chr", "", chromosome, ignore.case = TRUE)
  ) %>%
  # E??le??me anahtar??m??z?? olu??turuyoruz: "1:123456" format??nda
  mutate(join_key = paste(chromosome, base_pair_location, sep = ":"))

print("Lipo data haz??rland??.")

# ---------------------------------------------------------------------------------- #
# 2 - GWAS DATA STANDARDIZASYONU (t2d_suzuki_gwas)
# ---------------------------------------------------------------------------------- #
# Kolonlar??n: Chromsome, Position, EffectAllele, Pval, Neff...

gwas_clean <- t2d_suzuki_gwas %>%
  # ??nce ??u typo'yu ve isimleri d??zeltelim, kafam??z kar????mas??n
  dplyr::rename(
    chromosome = Chromsome,  # "Chromsome" -> "chromosome" yap??yoruz
    position = Position,
    p_value_gwas = Pval,     # Lipo'daki p-value ile kar????mas??n diye ad??n?? de??i??tirdim
    n_eff = Neff
  ) %>%
  # Sadece i??imize yarayacak kolonlar?? se??elim (RAM tasarrufu)
  dplyr::select(chromosome, position, EffectAllele, NonEffectAllele, Beta, SE, p_value_gwas, n_eff) %>%
  mutate(
    chromosome = as.character(chromosome),
    position = as.character(position),
    chromosome = gsub("chr", "", chromosome, ignore.case = TRUE)
  ) %>%
  # Ayn?? anahtar?? burada da olu??turuyoruz
  mutate(join_key = paste(chromosome, position, sep = ":"))

print("GWAS data haz??rland??.")

# ---------------------------------------------------------------------------------- #
# 3 - B??RLE??T??RME (MERGE / JOIN)
# ---------------------------------------------------------------------------------- #

# Lipo dosyas??n?? baz al??p, GWAS bilgilerini yan??na ??ekiyoruz.
merged_data <- lipo_clean %>%
  left_join(gwas_clean, by = "join_key")

# join_key ve duplicate kolonlar?? temizleyelim (chromosome.x, chromosome.y olmas??n)
merged_data <- merged_data %>%
  dplyr::select(-join_key, -chromosome.y, -position) %>%
  dplyr::rename(chromosome = chromosome.x)

print(paste("E??le??me tamamland??. Toplam sat??r:", nrow(merged_data)))

# ---------------------------------------------------------------------------------- #
# 4 - BIOMART ??LE EKS??K rsID'leri ??EKME (CHUNK Y??NTEM??)
# ---------------------------------------------------------------------------------- #

# Ba??lant?? kural??m
ensembl37 <- useMart(biomart = "ENSEMBL_MART_SNP", 
                     host    = "https://grch37.ensembl.org",
                     path    = "/biomart/martservice", 
                     dataset = "hsapiens_snp")

# Koordinatlar?? haz??rlayal??m
coords <- merged_data %>%
  dplyr::select(chromosome, base_pair_location) %>%
  dplyr::distinct() 

print(paste("Toplam sorgulanacak varyant say??s??:", nrow(coords)))

# --- CHUNKING BA??LIYOR ---

# 500'erli paketlere b??lelim (Server'?? k??zd??rmamak i??in g??venli say??)
batch_size <- 500 
chunks <- split(coords, ceiling(seq_along(rownames(coords))/batch_size))

results_list <- list() # Sonu??lar?? burada biriktirece??iz

print("Sorgu d??ng??s?? ba??l??yor, biraz sab??r...")

# ??lerleme ??ubu??u (Progress bar) - ne kadar kald??????n?? g??r diye
pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)

for (i in 1:length(chunks)) {
  
  chunk <- chunks[[i]]
  
  # Hata olursa d??ng?? k??r??lmas??n diye 'tryCatch' kullan??yoruz
  try({
    res <- getBM(
      attributes = c('refsnp_id', 'chr_name', 'chrom_start'), 
      filters    = c('chr_name', 'start', 'end'), 
      values     = list(chunk$chromosome, chunk$base_pair_location, chunk$base_pair_location), 
      mart       = ensembl37
    )
    
    results_list[[i]] <- res
  })
  
  # Server'a nefes ald??r (1 saniye bekle)
  Sys.sleep(1)
  
  # ??lerleme ??ubu??unu g??ncelle
  setTxtProgressBar(pb, i)
}

close(pb)

# T??m par??alar?? tek bir tablo haline getir
rsid_results <- dplyr::bind_rows(results_list)

print("Sorgu bitti! Veriler birle??tiriliyor...")

# --- VER?? ????LEME DEVAM ---

# Gelen sonu??lar?? string yapal??m
rsid_results <- rsid_results %>%
  mutate(
    chrom_start = as.character(chrom_start),
    # Join i??in anahtar
    join_key_bm = paste(chr_name, chrom_start, sep = ":")
  ) %>%
  dplyr::select(join_key_bm, refsnp_id)

# ---------------------------------------------------------------------------------- #
# 5 - F??NAL TABLO
# ---------------------------------------------------------------------------------- #

# Merged dataya rsID'leri ekliyoruz
final_table <- merged_data %>%
  mutate(join_key_bm = paste(chromosome, base_pair_location, sep = ":")) %>%
  left_join(rsid_results, by = "join_key_bm") %>%
  dplyr::select(-join_key_bm) %>% # Ara de??i??keni sil
  # E??er birden fazla rsID d??nerse, duplicate sat??r olu??abilir.
  # Bunu ??nlemek i??in rsID'leri virg??lle birle??tirebiliriz (Opsiyonel):
  group_by(chromosome, base_pair_location) %>%
  mutate(refsnp_id = paste(unique(refsnp_id), collapse = ",")) %>%
  distinct() %>%
  ungroup()

# Kontrol
print(head(final_table))

# Kaydet
write.csv(final_table, "T2D_Lipo_Full_with_GWAS_and_rsIDs.csv", row.names = FALSE)
print("????lem tamam kanka, dosya kaydedildi.")