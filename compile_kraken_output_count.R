library(wilkoxmisc)
#Open Table
LOR001 <- read.tidy("LOR001C_S8_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR001) <- c("Species","Abundance")
LOR001$Sample <- "LOR001"

LOR002 <- read.tidy("LOR002C_S15_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR002) <- c("Species","Abundance")
LOR002$Sample <- "LOR002"

LOR003 <- read.tidy("LOR003C_S24_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR003) <- c("Species","Abundance")
LOR003$Sample <- "LOR003"

LOR005 <- read.tidy("LOR005C_S17_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR005) <- c("Species","Abundance")
LOR005$Sample <- "LOR005"

LOR006 <- read.tidy("LOR006C_S18_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR006) <- c("Species","Abundance")
LOR006$Sample <- "LOR006"

LOR007 <- read.tidy("LOR007C_S19_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR007) <- c("Species","Abundance")
LOR007$Sample <- "LOR007"


LOR008 <- read.tidy("LOR008C_S20_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR008) <- c("Species","Abundance")
LOR008$Sample <- "LOR008"

LOR009 <- read.tidy("LOR009C_S29_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR009) <- c("Species","Abundance")
LOR009$Sample <- "LOR009"

LOR010 <- read.tidy("LOR010C_S23_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR010) <- c("Species","Abundance")
LOR010$Sample <- "LOR010"

LOR011 <- read.tidy("LOR011C_S23_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR011) <- c("Species","Abundance")
LOR011$Sample <- "LOR011"

LOR012 <- read.tidy("LOR012C_S25_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR012) <- c("Species","Abundance")
LOR012$Sample <- "LOR012"

LOR013 <- read.tidy("LOR013C_S33_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR013) <- c("Species","Abundance")
LOR013$Sample <- "LOR013"

LOR014 <- read.tidy("LOR014C_S26_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR014) <- c("Species","Abundance")
LOR014$Sample <- "LOR014"

LOR015 <- read.tidy("LOR015C_S27_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR015) <- c("Species","Abundance")
LOR015$Sample <- "LOR015"

LOR016 <- read.tidy("LOR016C_S36_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR016) <- c("Species","Abundance")
LOR016$Sample <- "LOR016"

LOR017 <- read.tidy("LOR017C_S30_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR017) <- c("Species","Abundance")
LOR017$Sample <- "LOR017"

LOR018 <- read.tidy("LOR018C_S31_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR018) <- c("Species","Abundance")
LOR018$Sample <- "LOR018"

LOR019 <- read.tidy("LOR019C_S31_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR019) <- c("Species","Abundance")
LOR019$Sample <- "LOR019"

LOR020 <- read.tidy("LOR020C_S40_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR020) <- c("Species","Abundance")
LOR020$Sample <- "LOR020"

LOR021 <- read.tidy("LOR021C_S41_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR021) <- c("Species","Abundance")
LOR021$Sample <- "LOR021"

LOR022 <- read.tidy("LOR022C_S42_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR022) <- c("Species","Abundance")
LOR022$Sample <- "LOR022"

LOR023 <- read.tidy("LOR023C_S43_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR023) <- c("Species","Abundance")
LOR023$Sample <- "LOR023"

LOR024 <- read.tidy("LOR024C_S37_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR024) <- c("Species","Abundance")
LOR024$Sample <- "LOR024"

LOR025 <- read.tidy("LOR025C_S45_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR025) <- c("Species","Abundance")
LOR025$Sample <- "LOR025"

LOR026 <- read.tidy("LOR026C_S9_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR026) <- c("Species","Abundance")
LOR026$Sample <- "LOR026"

LOR027 <- read.tidy("LOR027C_S39_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR027) <- c("Species","Abundance")
LOR027$Sample <- "LOR027"

LOR028 <- read.tidy("LOR028C_S47_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR028) <- c("Species","Abundance")
LOR028$Sample <- "LOR028"

LOR029 <- read.tidy("LOR029C_S48_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR029) <- c("Species","Abundance")
LOR029$Sample <- "LOR029"

LOR030 <- read.tidy("LOR030C_S42_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR030) <- c("Species","Abundance")
LOR030$Sample <- "LOR030"

LOR033 <- read.tidy("LOR033C_S50_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR033) <- c("Species","Abundance")
LOR033$Sample <- "LOR033"

LOR034 <- read.tidy("LOR034C_S51_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR034) <- c("Species","Abundance")
LOR034$Sample <- "LOR034"

LOR036 <- read.tidy("LOR036C_S52_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR036) <- c("Species","Abundance")
LOR036$Sample <- "LOR036"

LOR037 <- read.tidy("LOR037C_S53_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR037) <- c("Species","Abundance")
LOR037$Sample <- "LOR037"

LOR042 <- read.tidy("LOR042C_S47_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR042) <- c("Species","Abundance")
LOR042$Sample <- "LOR042"

LOR046 <- read.tidy("LOR046C_S48_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR046) <- c("Species","Abundance")
LOR046$Sample <- "LOR046"

LOR051 <- read.tidy("LOR051C_S56_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR051) <- c("Species","Abundance")
LOR051$Sample <- "LOR051"

LOR054 <- read.tidy("LOR054C_S50_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR054) <- c("Species","Abundance")
LOR054$Sample <- "LOR054"

LOR055 <- read.tidy("LOR055C_S58_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR055) <- c("Species","Abundance")
LOR055$Sample <- "LOR055"

LOR058 <- read.tidy("LOR058C_S52_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR058) <- c("Species","Abundance")
LOR058$Sample <- "LOR058"

LOR060 <- read.tidy("LOR060C_S60_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR060) <- c("Species","Abundance")
LOR060$Sample <- "LOR060"

LOR061 <- read.tidy("LOR061C_S53_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR061) <- c("Species","Abundance")
LOR061$Sample <- "LOR061"

LOR062 <- read.tidy("LOR062C_S55_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR062) <- c("Species","Abundance")
LOR062$Sample <- "LOR062"

LOR063 <- read.tidy("LOR063C_S56_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR063) <- c("Species","Abundance")
LOR063$Sample <- "LOR063"

LOR065 <- read.tidy("LOR065C_S57_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR065) <- c("Species","Abundance")
LOR065$Sample <- "LOR065"

LOR070 <- read.tidy("LOR070C_S58_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR070) <- c("Species","Abundance")
LOR070$Sample <- "LOR070"

LOR072 <- read.tidy("LOR072C_S59_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR072) <- c("Species","Abundance")
LOR072$Sample <- "LOR072"

LOR074 <- read.tidy("LOR074C_S60_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR074) <- c("Species","Abundance")
LOR074$Sample <- "LOR074"

LOR076 <- read.tidy("LOR076C_S68_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR076) <- c("Species","Abundance")
LOR076$Sample <- "LOR076"

LOR078 <- read.tidy("LOR078C_S10_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR078) <- c("Species","Abundance")
LOR078$Sample <- "LOR078"

LOR081 <- read.tidy("LOR081C_S61_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR081) <- c("Species","Abundance")
LOR081$Sample <- "LOR081"

LOR082 <- read.tidy("LOR082C_S63_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR082) <- c("Species","Abundance")
LOR082$Sample <- "LOR082"

LOR083 <- read.tidy("LOR083C_S64_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR083) <- c("Species","Abundance")
LOR083$Sample <- "LOR083"

LOR084 <- read.tidy("LOR084C_S65_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR084) <- c("Species","Abundance")
LOR084$Sample <- "LOR084"

LOR085 <- read.tidy("LOR085C_S65_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR085) <- c("Species","Abundance")
LOR085$Sample <- "LOR085"

LOR086 <- read.tidy("LOR086C_S67_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR086) <- c("Species","Abundance")
LOR086$Sample <- "LOR086"

LOR090 <- read.tidy("LOR090C_S75_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR090) <- c("Species","Abundance")
LOR090$Sample <- "LOR090"

LOR091 <- read.tidy("LOR091C_S68_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR091) <- c("Species","Abundance")
LOR091$Sample <- "LOR091"

LOR098 <- read.tidy("LOR098C_S69_R2_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR098) <- c("Species","Abundance")
LOR098$Sample <- "LOR098"

LOR099 <- read.tidy("LOR099C_S71_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR099) <- c("Species","Abundance")
LOR099$Sample <- "LOR099"

LOR101 <- read.tidy("LOR101C_S79_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR101) <- c("Species","Abundance")
LOR101$Sample <- "LOR101"

LOR102 <- read.tidy("LOR102C_S80_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR102) <- c("Species","Abundance")
LOR102$Sample <- "LOR102"

LOR103 <- read.tidy("LOR103C_S81_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR103) <- c("Species","Abundance")
LOR103$Sample <- "LOR103"

LOR202 <- read.tidy("LOR202C_S82_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR202) <- c("Species","Abundance")
LOR202$Sample <- "LOR202"

LOR203 <- read.tidy("LOR203C_S83_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR203) <- c("Species","Abundance")
LOR203$Sample <- "LOR203"

LOR204 <- read.tidy("LOR204C_S77_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR204) <- c("Species","Abundance")
LOR204$Sample <- "LOR204"

LOR206 <- read.tidy("LOR206C_S78_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR206) <- c("Species","Abundance")
LOR206$Sample <- "LOR206"

LOR207 <- read.tidy("LOR207C_S79_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR207) <- c("Species","Abundance")
LOR207$Sample <- "LOR207"

LOR208 <- read.tidy("LOR208C_S80_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR208) <- c("Species","Abundance")
LOR208$Sample <- "LOR208"

LOR210 <- read.tidy("LOR210C_S88_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR210) <- c("Species","Abundance")
LOR210$Sample <- "LOR210"

LOR211 <- read.tidy("LOR211C_S82_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR211) <- c("Species","Abundance")
LOR211$Sample <- "LOR211"

LOR212 <- read.tidy("LOR212C_S90_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR212) <- c("Species","Abundance")
LOR212$Sample <- "LOR212"

LOR214 <- read.tidy("LOR214C_S91_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR214) <- c("Species","Abundance")
LOR214$Sample <- "LOR214"

LOR215 <- read.tidy("LOR215C_S11_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR215) <- c("Species","Abundance")
LOR215$Sample <- "LOR215"

LOR216 <- read.tidy("LOR216C_S92_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR216) <- c("Species","Abundance")
LOR216$Sample <- "LOR216"

LOR217 <- read.tidy("LOR217C_S93_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR217) <- c("Species","Abundance")
LOR217$Sample <- "LOR217"

LOR218 <- read.tidy("LOR218C_S94_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR218) <- c("Species","Abundance")
LOR218$Sample <- "LOR218"

LOR219 <- read.tidy("LOR219C_S95_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR219) <- c("Species","Abundance")
LOR219$Sample <- "LOR219"

LOR220 <- read.tidy("LOR220C_S89_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR220) <- c("Species","Abundance")
LOR220$Sample <- "LOR220"

LOR221 <- read.tidy("LOR221C_S97_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR221) <- c("Species","Abundance")
LOR221$Sample <- "LOR221"

LOR222 <- read.tidy("LOR222C_S98_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR222) <- c("Species","Abundance")
LOR222$Sample <- "LOR222"

LOR223 <- read.tidy("LOR223C_S99_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR223) <- c("Species","Abundance")
LOR223$Sample <- "LOR223"

LOR224 <- read.tidy("LOR224C_S93_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR224) <- c("Species","Abundance")
LOR224$Sample <- "LOR224"

LOR225 <- read.tidy("LOR225C_S101_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR225) <- c("Species","Abundance")
LOR225$Sample <- "LOR225"

LOR226 <- read.tidy("LOR226C_S95_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR226) <- c("Species","Abundance")
LOR226$Sample <- "LOR226"

LOR227 <- read.tidy("LOR227C_S96_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR227) <- c("Species","Abundance")
LOR227$Sample <- "LOR227"

LOR228 <- read.tidy("LOR228C_S104_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR228) <- c("Species","Abundance")
LOR228$Sample <- "LOR228"

LOR229 <- read.tidy("LOR229C_S105_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR229) <- c("Species","Abundance")
LOR229$Sample <- "LOR229"

LOR230 <- read.tidy("LOR230C_S106_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR230) <- c("Species","Abundance")
LOR230$Sample <- "LOR230"

LOR231 <- read.tidy("LOR231C_S99_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR231) <- c("Species","Abundance")
LOR231$Sample <- "LOR231"

LOR233 <- read.tidy("LOR233C_S100_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR233) <- c("Species","Abundance")
LOR233$Sample <- "LOR233"

LOR234 <- read.tidy("LOR234C_S101_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR234) <- c("Species","Abundance")
LOR234$Sample <- "LOR234"

LOR235 <- read.tidy("LOR235C_S103_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR235) <- c("Species","Abundance")
LOR235$Sample <- "LOR235"

LOR237 <- read.tidy("LOR237C_S104_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR237) <- c("Species","Abundance")
LOR237$Sample <- "LOR237"

LOR238 <- read.tidy("LOR238C_S105_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR238) <- c("Species","Abundance")
LOR238$Sample <- "LOR238"

LOR239 <- read.tidy("LOR239C_S106_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR239) <- c("Species","Abundance")
LOR239$Sample <- "LOR239"

LOR240 <- read.tidy("LOR240C_S107_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR240) <- c("Species","Abundance")
LOR240$Sample <- "LOR240"

LOR241 <- read.tidy("LOR241C_S12_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR241) <- c("Species","Abundance")
LOR241$Sample <- "LOR241"

LOR242 <- read.tidy("LOR242C_S108_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR242) <- c("Species","Abundance")
LOR242$Sample <- "LOR242"

LOR243 <- read.tidy("LOR243C_S109_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR243) <- c("Species","Abundance")
LOR243$Sample <- "LOR243"

LOR246 <- read.tidy("LOR246C_S109_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR246) <- c("Species","Abundance")
LOR246$Sample <- "LOR246"

LOR248 <- read.tidy("LOR248C_S110_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR248) <- c("Species","Abundance")
LOR248$Sample <- "LOR248"

LOR249 <- read.tidy("LOR249C_S111_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR249) <- c("Species","Abundance")
LOR249$Sample <- "LOR249"

LOR254 <- read.tidy("LOR254C_S112_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR254) <- c("Species","Abundance")
LOR254$Sample <- "LOR254"

LOR258 <- read.tidy("LOR258C_S113_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR258) <- c("Species","Abundance")
LOR258$Sample <- "LOR258"

LOR259 <- read.tidy("LOR259C_S115_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR259) <- c("Species","Abundance")
LOR259$Sample <- "LOR259"

LOR260 <- read.tidy("LOR260C_S115_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR260) <- c("Species","Abundance")
LOR260$Sample <- "LOR260"

LOR262 <- read.tidy("LOR262C_S117_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR262) <- c("Species","Abundance")
LOR262$Sample <- "LOR262"

LOR264 <- read.tidy("LOR264C_S118_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR264) <- c("Species","Abundance")
LOR264$Sample <- "LOR264"

LOR269 <- read.tidy("LOR269C_S119_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR269) <- c("Species","Abundance")
LOR269$Sample <- "LOR269"

LOR270 <- read.tidy("LOR270C_S119_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR270) <- c("Species","Abundance")
LOR270$Sample <- "LOR270"

LOR272 <- read.tidy("LOR272C_S121_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR272) <- c("Species","Abundance")
LOR272$Sample <- "LOR272"

LOR274 <- read.tidy("LOR274C_S121_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR274) <- c("Species","Abundance")
LOR274$Sample <- "LOR274"

LOR275 <- read.tidy("LOR275C_S122_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR275) <- c("Species","Abundance")
LOR275$Sample <- "LOR275"

LOR279 <- read.tidy("LOR279C_S123_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR279) <- c("Species","Abundance")
LOR279$Sample <- "LOR279"

LOR281 <- read.tidy("LOR281C_S125_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR281) <- c("Species","Abundance")
LOR281$Sample <- "LOR281"

LOR282 <- read.tidy("LOR282C_S126_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR282) <- c("Species","Abundance")
LOR282$Sample <- "LOR282"

LOR283 <- read.tidy("LOR283C_S126_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR283) <- c("Species","Abundance")
LOR283$Sample <- "LOR283"

LOR284 <- read.tidy("LOR284C_S127_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR284) <- c("Species","Abundance")
LOR284$Sample <- "LOR284"

LOR287 <- read.tidy("LOR287C_S129_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR287) <- c("Species","Abundance")
LOR287$Sample <- "LOR287"

LOR288 <- read.tidy("LOR288C_S129_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR288) <- c("Species","Abundance")
LOR288$Sample <- "LOR288"

LOR289 <- read.tidy("LOR289C_S13_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR289) <- c("Species","Abundance")
LOR289$Sample <- "LOR289"

LOR294 <- read.tidy("LOR294C_S130_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR294) <- c("Species","Abundance")
LOR294$Sample <- "LOR294"

LOR295 <- read.tidy("LOR295C_S132_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR295) <- c("Species","Abundance")
LOR295$Sample <- "LOR295"

LOR300 <- read.tidy("LOR300C_S132_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR300) <- c("Species","Abundance")
LOR300$Sample <- "LOR300"

LOR301 <- read.tidy("LOR301C_S133_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR301) <- c("Species","Abundance")
LOR301$Sample <- "LOR301"

LOR303 <- read.tidy("LOR303C_S134_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(LOR303) <- c("Species","Abundance")
LOR303$Sample <- "LOR303"

NC1 <- read.tidy("NC1_S14_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(NC1) <- c("Species","Abundance")
NC1$Sample <- "NC1"

NC2 <- read.tidy("NC2_S135_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(NC2) <- c("Species","Abundance")
NC2$Sample <- "NC2"

NC3 <- read.tidy("NC3_S136_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(NC3) <- c("Species","Abundance")
NC3$Sample <- "NC3"

NC4 <- read.tidy("NC4_S137_R1_001_kneaddata_paired_1.fastq_rep.txt_bracken.txt")[c("name","new_est_reads")]
colnames(NC4) <- c("Species","Abundance")
NC4$Sample <- "NC4"

Combined <-

Combined <- rbind(LOR001,LOR002,LOR003,LOR005,LOR006,LOR007,LOR008,LOR009,LOR010,LOR011,LOR012,LOR013,LOR014,LOR015,LOR016,LOR017,LOR018,LOR019,LOR020,
                    LOR021,LOR022,LOR023,LOR024,LOR025,LOR026,LOR027,LOR028,LOR029,LOR030,LOR033,LOR034,LOR036,LOR037,LOR042,LOR046,
                    LOR051,LOR054,LOR055,LOR058,LOR060,LOR061,LOR062,LOR063,LOR065,LOR070,LOR072,LOR074,LOR076,LOR078,LOR081,LOR082,LOR083,LOR084,
                    LOR085,LOR086,LOR090,LOR091,LOR098,LOR099,LOR101,LOR102,LOR103,LOR202,LOR203,LOR204,LOR206,LOR207,LOR208,LOR210,LOR211,LOR212,
                    LOR214,LOR215,LOR216,LOR217,LOR218,LOR219,LOR220,LOR221,LOR222,LOR223,LOR224,LOR225,LOR226,LOR227,LOR228,LOR229,LOR230,LOR231,
                    LOR233,LOR234,LOR235,LOR237,LOR238,LOR239,LOR240,LOR241,LOR242,LOR243,LOR246,LOR248,LOR249,LOR254,LOR258,LOR259,LOR260,LOR262,
                    LOR264,LOR269,LOR270,LOR272,LOR274,LOR275,LOR279,LOR281,LOR282,LOR283,LOR284,LOR287,LOR288,LOR289,LOR294,LOR295,LOR300,LOR301,LOR303,
                    NC1,NC2,NC3,NC4)

write.tidy(Combined,"bracken2_all_samples_long_count.txt")

#Convert from long to wide
library(dplyr)
Wide <- dcast(Combined,Species~Sample, fun.agg = function(Combined) sum(!is.na(Combined)), value.var = "Abundance")


#Draw species-specific plots of abundance between Baoding and Dalian
Table <- read.tidy("bracken2_all_samples_wide_0_00001_adjusted_for_human_clean.txt")
Table <- melt(Table)
Meta <- read.tidy("../meta/meta_new.txt")[c("Sample","City")]
Merge <- merge(Table,Meta,by.x="variable",by.y="Sample",all=TRUE)

#Filter table by species
C_acnes <- Merge[which(Merge$Species == "Cutibacterium acnes") ,]
Plot <- ggplot(C_acnes,aes(x=City,y=value,fill=City))
Plot <- Plot + geom_boxplot() + theme_classic()
Plot <- Plot + xlab(paste0("City")) + ylab(paste0("Relative Abundance (%)")) + theme(legend.position="none")
ggsave("city_acnes.pdf",width=6,height=6,units="in")
#stat
wilcox.test(value~City,data=C_acnes)
#data:  value by City
#W = 5.231e+09, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

S_aureus <- Merge[which(Merge$Species == "Staphylococcus aureus") ,]
Plot <- ggplot(S_aureus,aes(x=City,y=value,fill=City))
Plot <- Plot + geom_boxplot(outlier.shape = NA)  + theme_classic()
Plot <- Plot + xlab(paste0("City")) + ylab(paste0("Relative Abundance (%)")) + theme(legend.position="none") + scale_y_continuous(limits = c(0, 20)) 
ggsave("city_Saureus.pdf",width=6,height=6,units="in")
#stat
wilcox.test(value~City,data=S_aureus)
#data:  value by City
#W = 5.231e+09, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

M_luteus <- Merge[which(Merge$Species == "Micrococcus luteus") ,]
Plot <- ggplot(M_luteus,aes(x=City,y=value,fill=City))
Plot <- Plot + geom_boxplot(outlier.shape = NA)  + theme_classic()
Plot <- Plot + xlab(paste0("City")) + ylab(paste0("Relative Abundance (%)")) + theme(legend.position="none") + scale_y_continuous(limits = c(0, 10)) 
ggsave("city_Mluteus.pdf",width=6,height=6,units="in")
#stat
wilcox.test(value~City,data=M_luteus)

P_yeei <- Merge[which(Merge$Species == "Paracoccus yeei") ,]
Plot <- ggplot(P_yeei,aes(x=City,y=value,fill=City))
Plot <- Plot + geom_boxplot(outlier.shape = NA)  + theme_classic()
Plot <- Plot + xlab(paste0("City")) + ylab(paste0("Relative Abundance (%)")) + theme(legend.position="none") + scale_y_continuous(limits = c(0, 2.5)) 
ggsave("city_Pyeei.pdf",width=6,height=6,units="in")
#stat
wilcox.test(value~City,data=P_yeei)

S_epidermidis <- Merge[which(Merge$Species == "Staphylococcus epidermidis") ,]
Plot <- ggplot(S_epidermidis,aes(x=City,y=value,fill=City))
Plot <- Plot + geom_boxplot(outlier.shape = NA)  + theme_classic()
Plot <- Plot + xlab(paste0("City")) + ylab(paste0("Relative Abundance (%)")) + theme(legend.position="none") 
ggsave("city_Sepidermidis.pdf",width=6,height=6,units="in")
#stat
wilcox.test(value~City,data=S_epidermidis)

#Visualize C. acnes and phage relationships by city
Table <- read.tidy("cacnes_phage_table.txt")
Plot <- ggplot(Table,aes(x=C_acnes_per,y=Phage_per,colour=City))
Plot <- Plot + geom_point() + geom_smooth(size=0.5, method=lm, se=FALSE, colour="black", formula = y ~ x) + facet_wrap(~City,scales="free")
Plot <- Plot + theme_classic() + xlab(paste0("C. acnes Relative Abundance (%)")) + ylab(paste0("C. acnes Phage Relative Abundance (%)"))
ggsave("C_acnes_phage_relation.pdf",width=6,height=9,units="in")
Bao <- Table[which(Table$City == "Baoding") ,]
cor.test(Bao$C_acnes_per,Bao$Phage_per,data=Bao,method="pearson")
#data:  Bao$C_acnes_per and Bao$Phage_per
#t = 2.208, df = 58, p-value = 0.03121
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.02639624 0.49722241
#sample estimates:
#  cor 
#0.2784546 

Dal <- Table[which(Table$City == "Dalian") ,]
cor.test(Dal$C_acnes_per,Dal$Phage_per,data=Dal,method="pearson")

#Draw taxonomy plot for top organisms
library(devtools)
library(wilkoxmisc)
library(reshape2)
library(ggplot2)

#Open taxonomy OTU table
Species <- read.tidy("bracken2_all_samples_wide_0_00001_adjusted_for_human_clean_no_con.txt")
Species <- melt(Species)
names(Species)[2] <- "Sample"
names(Species)[3] <- "Abundance"
#Open Metadata table
Meta <- read.tidy("../meta/meta_new.txt")[c("Sample","City")]

#Tabulate read counts by genus
Species <- ddply(Species, .(Sample, Species), summarise, Abundance = sum(Abundance))
#Convert count to relativeabundance and add column
Species <- ddply(Species, .(Sample), mutate, RelativeAbundance = (Abundance * 100) / sum(Abundance))

#collapse taxa table to only 5 or 8 top phyla, genus, family, etc (require reshape2).
SpeciesTable <- collapse.taxon.table(Species, n = 15, Rank = "Species")

#Merge relative abundance table and metatable together
SpeciesTable <- merge(SpeciesTable, Meta, by = "Sample", all.x = TRUE)

write.tidy(OTUTable, "Top15_taxonomy_species.txt")

#Construct histogram of prevalence count
Table <- read.tidy("prevalence_count_histogram_plot.txt")

Plot <- ggplot(Table, aes(x=Prevalence,y=Count))
Plot <- Plot + geom_bar(stat="identity") + theme_classic() + xlab(paste0("Number of Samples")) +ylab(paste0("Number of Species")) + scale_y_continuous(expand=c(0,0))
ggsave("histogram_prevalence.pdf",width=6,height=6,units="in")