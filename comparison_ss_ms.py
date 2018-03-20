CFViewBackward(1210)
FileOpenProject('C:/Users/msais/Box Sync/Thesis Work/Baseline/baseline2/baseline2_4kgs_SA_mav/baseline2_4kgs_SA_mav.run')
ViewActivate('baseline2_4kgs_SA_mav.run:1')
SetNumecaLogo(0, 0)
ViewOpenSTM(0.238415,0.338415,0.42454,0.52454)
Limits(-1.00305,-0.00304878,-0.97546,1.02454)
ViewActivate('baseline2_4kgs_SA_mav.run:1')
ViewClose()
ViewActivate('baseline2_4kgs_SA_mav.run:2')
GmtRepetitionToggle()
ViewAxes()

FileOpenProject('C:/Users/msais/Box Sync/Thesis Work/Multi-Stage_data/DiffuserConstArea/WorkSplitRotor1=35/Stator25deg/MSD_Beta/MSD_Beta_4kgs_111_mav/MSD_Beta_4kgs_111_mav.run')
SetNumecaLogo(0, 0)
ViewTile()
ViewActivate('MSD_Beta_4kgs_111_mav.run:3')
ViewOpenSTM(0.238415,0.338415,0.42454,0.52454)
Limits(-1.00305,-0.00304878,-0.97546,1.02454)
ViewActivate('MSD_Beta_4kgs_111_mav.run:3')
ViewClose()
ViewActivate('MSD_Beta_4kgs_111_mav.run:4')
GmtRepetitionToggle()
ViewAxes()
file_dir = 'C:/Users/msais/Box Sync/Thesis Work/Report/Comparing SS_MS contour pics/'

Qnt_vec = ['Vt']
span_vec = [0.25, 0.50, 0.75, 0.90]

EntropyRange = [-300, 10]

cntr=0
for Qnt in Qnt_vec:
    for span in span_vec:
        ViewActivate('MSD_Beta_4kgs_111_mav.run:4')
        GmtToggleBoundary()
        CutPlaneSave(span,0,0,1,0,0,2,0,str(span))
        GmtToggleBoundary()
        SelectFromProject(str(span))
        QntFieldScalar(Qnt)
        SclContourStrip()
        ColormapStripesOnly()
        ColormapTicksNumberTextType(10,16,2,0,1,0,0,1,0,0,0,0)
        ColormapLabelTextType(10,18,2,2,1,0,0,1,0,0,0,0)
        RprRangeIn(EntropyRange[0], EntropyRange[1])
        ViewZoomAll(1)

        ViewActivate('baseline2_4kgs_SA_mav.run:2')
        GmtToggleBoundary()
        CutPlaneSave(span,0,0,1,0,0,2,0,str(span))
        GmtToggleBoundary()
        SelectFromProject(str(span))
        QntFieldScalar(Qnt)
        SclContourStrip()
        ColormapStripesOnly()
        ColormapTicksNumberTextType(10,16,2,0,1,0,0,1,0,0,0,0)
        ColormapLabelTextType(10,18,2,2,1,0,0,1,0,0,0,0)
        ViewZoomAll(1)
        ViewTile()
        MatchRanges('MSD_Beta_4kgs_111_mav.run:4')
        Print(8, 1, 1, 1, 100, 1920, 1080, 0, file_dir + Qnt + '_B2b_' +
                str(span) + '.png', '', 1, 1, 1)
        DeleteAll()
        ViewActivate('MSD_Beta_4kgs_111_mav.run:4')
        DeleteAll()
        cntr = cntr + 1
