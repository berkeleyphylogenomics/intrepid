Ext.onReady(function(){
    
    var inform = new Ext.form.FormPanel({
       
        standardSubmit: true,
        frame: true,
        title: 'submit FASTA-formatted protein sequence',
        width: 600,
        defaultType: 'textarea',
        
        items: [{
            width: 540,
            height: 150,
            fieldLabel: 'sequence',
            name: 'sequence',
            id: 'sequence',
            allowBlank: true,
            blankText: 'A FASTA-formatted protein sequence is required.',
            style: 'font-size: small; font-family: monospace;'
        },
        new Ext.form.TextField({
            name: 'email',
            id: 'email',
            fieldLabel: 'email address (optional)',
            width: 540,
            allowBlank: true,
            vtype: 'email',
            vtypeText: 'email address must be of form: "user@domain.edu"'
        })],
            
        buttons: [{
            text: 'Submit',
            handler: function() {
                inform.getForm().getEl().dom.action = '../../cgi-bin/intrepid/run_intrepid.py';
                inform.getForm().getEl().dom.method = 'POST';
                
                var emailbox = inform.findById('email');

                if(emailbox.isValid()) {
                    inform.getForm().submit();
                }
                else {
                    Ext.Msg.alert('Invalid Email Address', 'email address (if provided) must valid, like: "user@domain.edu"');
                }
            }
        },
        {
            text: 'Clear',
            handler: function() {
                var seqarea = inform.findById('sequence');
                seqarea.setValue('');
            }
        },
        {
            text: 'Load Example 1',
            handler: function() {
                var seqarea = inform.findById('sequence');
                seqarea.setValue('>pdb|2JEK|A Chain A, Crystal Structure Of The Conserved Hypothetical Protein Rv1873 From Mycobacterium Tuberculosis At 1.38 A\nMKSASDPFDLKRFVYAQAPVYRSVVEELRAGRKRGHWMWFVFPQLRGLGSSPLAVRYGIS\nSLEEAQAYLQHDLLGPRLHECTGLVNQVQGRSIEEIFGPPDDLKLCSSMTLFARATDANQ\nDFVALLAKYYGGGEDRRTVALLAVT');
            }
        },
        {
            text: 'Load Example 2',
            handler: function() {
                var seqarea = inform.findById('sequence');
                seqarea.setValue('>Sporozoite stage TSP1 domain protein, S21\nMLMKISRYFFLLYLIKAHLDFFLRYRTGFIRSRLETYIGNSDVRYNKSFINNRLLNEHAH\nCDAWSEWSACSKTCDYGIKIRVKISTDQTKSKACSNITESTICHEHICPRTFEEAEETYL\nHNKEKEKKKKFRTTYILIFTIFSVVNIVVLLICVILSIKKKII');
            }
        }]
    });

    inform.render('seqinputform');

});
