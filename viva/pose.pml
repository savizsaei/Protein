load /Users/saviz/Documents/Saviz/SelfStudy/Protein/viva/receptor.pdb, receptor
load /Users/saviz/Documents/Saviz/SelfStudy/Protein/viva/docked_aspirin.pdb, lig
hide everything
show cartoon, receptor
color slate, receptor
show sticks, lig
color yellow, lig
bg_color white
set ray_opaque_background, off
zoom lig, 8
ray 1600,1200
png /Users/saviz/Documents/Saviz/SelfStudy/Protein/viva/pose.png, dpi=200
