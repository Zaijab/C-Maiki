function Heatmap(M)
xvalues = {'Species 1','Species 2','Species 3','Species 4','Species 5','Species 6','Species 7','Species 8','Species 9',...
    'Species 10','Species 11','Species 12','Species 13','Species 14','Species 15','Species 16','Species 17','Species 18','Species 19',...
    'Species 20','Species 21','Species 22','Species 23','Species 24','Species 25','Species 26','Species 27','Species 28','Species 29',...
    'Species 30','Species 31','Species 32','Species 33','Species 34','Species 35','Species 36','Species 37','Species 38','Species 39',...
    'Species 40','Species 41','Species 42','Species 43','Species 44','Species 45','Species 46','Species 47','Species 48','Species 49',...
    'Species 50','Species 51','Species 52','Species 53','Species 54','Species 55','Species 56','Species 57','Species 58','Species 59',...
    'Species 60','Species 61','Species 62','Species 63','Species 64','Species 65','Species 66','Species 67','Species 68','Species 69',...
    'Species 70','Species 71','Species 72','Species 73','Species 74','Species 75','Species 76','Species 77','Species 78','Species 79',...
    'Species 80','Species 81','Species 82','Species 83','Species 84','Species 85','Species 86','Species 87','Species 88','Species 89',...
    'Species 90','Species 91','Species 92','Species 93','Species 94','Species 95','Species 96','Species 97','Species 98','Species 99',...
    'Species 100','Species 101','Species 102','Species 103','Species 104','Species 105','Species 106','Species 107','Species 108','Species 109',...
    'Species 110','Species 111','Species 112','Species 113','Species 114','Species 115','Species 116','Species 117','Species 118','Species 119',...
    'Species 120','Species 121','Species 122','Species 123','Species 124','Species 125','Species 126','Species 127','Species 128','Species 129',...
    'Species 130','Species 131','Species 132','Species 133','Species 134','Species 135','Species 136','Species 137','Species 138','Species 139',...
    'Species 140','Species 141','Species 142','Species 143','Species 144','Species 145','Species 146','Species 147','Species 148','Species 149',...
    'Species 150','Species 151','Species 152','Species 153','Species 154','Species 155','Species 156','Species 157','Species 158','Species 159',...
    'Species 160','Species 161','Species 162','Species 163','Species 164','Species 165','Species 166','Species 167','Species 168','Species 169',...
    'Species 170','Species 171','Species 172','Species 173','Species 174','Species 175','Species 176','Species 177','Species 178','Species 179',...
    'Species 180','Species 181','Species 182','Species 183','Species 184','Species 185','Species 186','Species 187','Species 188','Species 189',...
    'Species 190','Species 191','Species 192','Species 193','Species 194','Species 195','Species 196','Species 197','Species 198','Species 199',...
    'Species 200','Species 201','Species 202','Species 203','Species 204','Species 205','Species 206','Species 207','Species 208','Species 209',...
    'Species 210','Species 211','Species 212','Species 213','Species 214','Species 215','Species 216','Species 217','Species 218','Species 219',...
    'Species 220','Species 221','Species 222','Species 223','Species 224','Species 225','Species 226','Species 227','Species 228','Species 229',...
    'Species 230','Species 231','Species 232','Species 233','Species 234','Species 235','Species 236','Species 237','Species 238','Species 239',...
    'Species 240','Species 241','Species 242','Species 243','Species 244','Species 245','Species 246','Species 247','Species 248','Species 249',...
    'Species 250','Species 251','Species 252','Species 253','Species 254','Species 255','Species 256','Species 257','Species 258','Species 259',...
    'Species 260','Species 261','Species 262','Species 263','Species 264','Species 265','Species 266','Species 267','Species 268','Species 269',...
    'Species 270','Species 271','Species 272','Species 273','Species 274','Species 275','Species 276','Species 277','Species 278','Species 279',...
    'Species 280','Species 281','Species 282','Species 283','Species 284','Species 285','Species 286','Species 287','Species 288','Species 289',...
    'Species 290','Species 291','Species 292','Species 293'};
yvalues={'Tree 1','Tree 2','Tree 3','Tree 4','Tree 5','Soil 1','Tree 6'...
    'Tree 7', 'Tree 8', 'Tree 9', 'Water 1', 'Soil 2', 'Water 2',...
    'Tree 10','Soil 3', 'Road 1', 'Road 2', 'Tree 11', 'Tree 12',...
    'Tree 13'};
h=heatmap(xvalues,yvalues,M);
h.XLabel = 'Species';
h.YLabel = 'Site';