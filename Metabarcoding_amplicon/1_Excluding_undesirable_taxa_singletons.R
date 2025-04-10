###

# Set colors for plotting
phylum_colors <- c(
  "gray",'black', "#DA5724", "#5F7FC7","#508578", "#CD9BCD", "orange",
  "#5F7FC7","#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#D1A33D", "#8A7C64", "#599861","#5E738F"
)

my_color_collection <- c(
  "#CBD588", "#5F7FC7", "orange", "#AD6F3B", "#673770",
  "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F",
  "#D1A33D", "#8A7C64", "#599861","#616163", "#FFCDB2",
  "#6D9F71", "#242F40",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")

my_color_Class <- c(
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "orange", "#5F7FC7", "#CBD588", "#AD6F3B", "#673770",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',

  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")



my_color_Order <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

my_color_Family <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)
my_color_Family2 <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

# up to 300
my_color_OTU <- c(
  "gray",'black', "#5F7FC7", "orange",  "#AD6F3B",
  "#673770","#D14285", "#652926", "#C84248",  "#8569D5",
  "#5E738F","#D1A33D", "#8A7C64", "#599861","#616163",
  "#FFCDB2", "#242F40", "#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482',
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy',
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2', '#FFF275', 'springgreen', '#FF3C38', '#A23E48',
  '#000000', '#CF5C36', '#EEE5E9', '#7C7C7C', '#EFC88B',

  '#2E5266', '#6E8898', '#9FB1BC', '#D3D0CB', '#E2C044',
  '#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921',
  "#CD9BCD", "#508578", "#CBD588","#CBD588", "#5F7FC7",
  "orange",   "#AD6F3B", "#673770","#D14285", "#652926",
  "#C84248",  "#8569D5", "#5E738F","#D1A33D", "#8A7C64",
  "#599861","#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", "#DEC4A1",
  "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', '#95C623',
  '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', '#A5668B',
  '#69306D', '#0E103D', '#1A535C', '#4ECDC4', '#F7FFF7',
  '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange")

my_color_otu2 <- c(
  "gray",'black', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB', 'pink',
  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "orange")


my_color_gen <- c(
  "white",'yellow', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#CD9BCD', '#6699CC', 'pink',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB',

  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357')


### For bacteria ####
tax_table(phy)
#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
phy.cp <- subset_taxa(phy, Order == "o__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 30 ASVs of CP
vec.cp #30
# (2) MT
phy.mt1 <- subset_taxa(phy, Family == "f__Mitochondria") #16 (latest DB)
#phy.mt <- subset_taxa(phy, Order == "o__Rickettsiales") #17 ASVs (latest DB)
vec.mt <- rownames(otu_table(phy.mt1))
tax_table(phy.mt1)
length(rownames(otu_table(phy.mt1))) ## 14 ASVs of MT

# (3) Unassigned and Archaea
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom %in% c("Unassigned","d__Archaea"))
tax_table(phy.un)
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 78 ASVs (Latest DB)

### exclude those vectors
phy  #7967 taxa and 54 samples

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

bac.clean <- pop_taxa(phy, c(vec.cp, vec.mt, vec.un))

##after clean-up
bac.clean ##7844
sum(otu_table(bac.clean)) #1497473
sort(colSums(otu_table(bac.clean))) #minimum 18472 reads, maximum 56533 reads

# checking procedure of whether the MT and CP otus are cleaned
#taxa_names(subset_taxa(phy, Order == "o__Rickettsiales"))  # before 10
taxa_names(subset_taxa(phy , Family=="f__Mitochondria")) # after 0
taxa_names(subset_taxa(phy, Order=="o__Chloroplast")) # before 305

#taxa_names(subset_taxa(bac.clean , Order == "o__Rickettsiales")) # after 0
taxa_names(subset_taxa(bac.clean , Family=="f__Mitochondria")) # after 0
taxa_names(subset_taxa(bac.clean , Order == "o__Chloroplast")) # after 0



#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt()
# phy.test4

tax_table(bac.clean)[,colnames(tax_table(bac.clean))] <- gsub(tax_table(bac.clean)[,colnames(tax_table(bac.clean))],pattern="[a-z]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(bac.clean)

### filter otu with total count of 20? (in all samples)
### later we can implement
bac.clean.asv <- otu_table(bac.clean)
head(bac.clean.asv)
df.clean.asv <- data.frame(bac.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(bac.clean)

## Remove reads with over 350 bp
library(seqinr)
bac.seq <- read.fasta(file = "./Bacteria/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
hist(getLength(bac.seq))

asv_over_350bp <- attr(bac.seq[which(getLength(bac.seq)>350)], "names")
asv_less_250bp <- attr(bac.seq[which(getLength(bac.seq)<250)], "names")


bac.clean
sum(otu_table(bac.clean)) #1512840
bac.clean.ss <- pop_taxa(bac.clean,asv_over_350bp)
bac.clean.ss<- pop_taxa(bac.clean.ss,asv_less_250bp)
bac.clean.ss  ## 7843 ASVs (Latest DB)
sum(otu_table(bac.clean.ss)) # 1497470 (Latest DB)
sort(colSums(otu_table(bac.clean.ss))) #minimum 18472 reads, maximum 56533 reads


## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
bac.list$number <- paste0('B',1:dim(bac.list)[1])

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

bac.list
##### For fungi ####
# (3) Unassigned
unique(tax_table(fun)[,'Kingdom'])
tax_table(fun)[,'Kingdom']
unique(tax_table(fun)[,'Kingdom'])
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Eukaryota_kgd_Incertae_sedis","k__Viridiplantae","k__Rhizaria","k__unidentified", "k__Protista"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ## 121 ASVs



### exclude those vectors
fun  # 54 samples,  2082 taxa

# sample_names(fun)
# sample_variables(fun)
# 
# ### pop taxa application
# ## get rid of CP and MT otus
# ### pop_taxa function
# pop_taxa = function(physeq, badTaxa){
#   allTaxa = taxa_names(physeq)
#   myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
#   return(prune_taxa(myTaxa, physeq))
# }
# ### Clean it up!!!
fun.clean <- pop_taxa(fun, c(vec.un))
#### We will also remove the "D_3__" patterns for cleaner labels
tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean    ## 1961 taxa
str(fun.clean)
otu_table(fun.clean)

fun.clean 


# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement
fun.clean.asv <- otu_table(fun.clean)
head(fun.clean.asv)
df.clean.asv <- data.frame(fun.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(fun.clean)

##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "./Fungi/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
hist(getLength(fun.seq))

asv_less_than_200bp <- attr(fun.seq[which(getLength(fun.seq)<200)], "names")
fun.clean.ss <- pop_taxa(fun.clean,asv_less_than_200bp)

fun.clean.ss ##  1944
colSums(otu_table(fun.clean.ss)) 


## Designating OTU id
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id

fun.list

otu.list <- rbind(bac.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id



#### Split DNA and RNA data for the bacterial data
bac.clean.ss.dna <- subset_samples(bac.clean.ss, NucleicAcids == "DNA")
bac.clean.ss.dna <- phyloseq::filter_taxa(bac.clean.ss.dna, function(x) sum(x) != 0, TRUE) #7502 taxa

### Compare DNA with RNA samples
bac.clean.ss.rna <- subset_samples(bac.clean.ss, sampleID %in% c("cd1","cd2","cd3","cd4", "j1", "k1","k10","t1"))
bac.clean.ss.rna <- phyloseq::filter_taxa(bac.clean.ss.rna, function(x) sum(x) != 0, TRUE) #1444 taxa




##### Summary
otu.tab.bac<-otu_table(bac.clean.ss.dna)
otu.tab.bac[otu.tab.bac>0] <- 1
head(otu.tab.bac)
num.read<-data.frame(colSums(otu_table(bac.clean.ss.dna)))
names(num.read)[1] <- "Read"
num.asv<-data.frame(colSums(otu.tab.bac))
names(num.asv)[1] <- "Number_of_ASV"

summary.asv.all <- cbind(num.asv, num.read)

map <- sample_data(bac.clean.ss.dna)
for (i in rownames(map)){
  map$Read[which(rownames(map)==i)] <- summary.asv.all$Read[which(rownames(summary.asv.all) == i)]
  map$Number_of_ASV[which(rownames(map)==i)] <- summary.asv.all$Number_of_ASV[which(rownames(summary.asv.all) == i)]
}

head(map)

write.csv(map, "230724_number of bacterial ASVs and reads in all samples.csv")


otu.tab.fun<-otu_table(fun.clean.ss.f)
otu.tab.fun[otu.tab.fun>0] <- 1
head(otu.tab.fun)
num.read<-data.frame(colSums(otu_table(fun.clean.ss.f)))
names(num.read)[1] <- "Read"
num.asv<-data.frame(colSums(otu.tab.fun))
names(num.asv)[1] <- "Number_of_ASV"

summary.asv.all <- cbind(num.asv, num.read)

map <-sample_data(fun.clean.ss.f)
for (i in rownames(map)){
  map$Read[which(rownames(map)==i)] <- summary.asv.all$Read[which(rownames(summary.asv.all) == i)]
  map$Number_of_ASV[which(rownames(map)==i)] <- summary.asv.all$Number_of_ASV[which(rownames(summary.asv.all) == i)]
}

head(map)

write.csv(map, "230724_number of fungal ASVs and reads in all samples.csv")


bac.clean.ss.dna
summarize_phyloseq(bac.clean.ss.dna)
summarize_phyloseq(fun.clean.ss)


save(bac.clean.ss.dna, file="~/Desktop/SNU/Conference & Seminar/KSPP/2024/PhytobiomeWorkshop/Test set 2/T_matsutkae_soil_bac_phyloseq.RData")
